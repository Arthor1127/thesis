from collections import OrderedDict

import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.qubits.squid_array import SquidArray

# ---------------------------------------------------------------------------
# This is the geometry drawn interactively in build_design.ipynb: a chain of
# N SQUID loops wired in series (a tunable superinductor - SquidArray) sitting
# next to a meandered readout resonator that runs between an open and a short
# termination. There is no transmon: the resonator's current threads flux
# through every loop, and an eigenmode EPR run reads how much of each mode's
# inductive energy sits in the array's junctions (see build_and_eigenmode.py
# and ../squid_resonator_epr/EPR_C_EXTRACTION.md for the method).
# ---------------------------------------------------------------------------

# ---- SQUID array: junction electrical parameters ----
N_SQUIDS = 23
L_J_NH = 15.0        # inductance of EVERY junction in the array (SquidArray.L_j)
CJ_J_FF = 1.5        # capacitance of every junction (SquidArray.C_j)
JJ_WIDTH_UM = 2.0    # transverse width of each junction bridge

# ---- SQUID array: single-loop geometry (um), as passed in the notebook ----
SQUID_GAP_UM = 25.0        # -> squid_gap        (loop height)
ARM_WIDTH_UM = 4.0         # -> plate1_width / segment_a_width / segment_b_width / segment_c_width
JOSEPHSON_GAP_UM = 3.0     # -> JJ_gap           (the junction sits in this gap)
SEGMENT_LENGTH_UM = 6.0    # -> segment_d_length (right contact length)
ARM_SEMI_LENGTH_UM = 15.0  # -> segment_a_length / segment_b_length
LINE_WIDTH_UM = 9.0        # -> segment_d_width  (right contact height)
GROUND_GAP_UM = 3.0        # -> ground_gap       (ground-plane pocket around the chain)

# ---- SQUID array: placement ----
ARRAY_POS_X_UM = -1290.0
ARRAY_POS_Y_UM = 0.0
ARRAY_ORIENTATION_DEG = 0.0
ARRAY_COMPONENT_NAME = 'array_1'

# SQDMetal's create_port_JosephsonJunction addresses junctions on a component
# positionally (by insertion order into design.qgeometry.tables['junction']),
# not by name. SquidArray emits, per loop i: jj_i_upper then jj_i_lower. So
# loop i -> junction_index 2*i (upper) and 2*i + 1 (lower); the array carries
# 2 * N_SQUIDS junction rows in total. Every gap is a real vacuum break, so
# every one MUST get a lumped inductive port - an un-ported gap is an open
# circuit that kills the loop's circulating current (EPR_C_EXTRACTION.md Sec 2).
JJ_PER_LOOP = 2

# ---- readout resonator ----
RES_LENGTH_MM = 3.74575

# ---- backward-compat aliases ----
# build_design.ipynb imports these names and is not to be edited. In this
# design the only inductive elements are the array junctions, so the old
# "probe" names alias the array junction values; the "transmon" names are
# unused (there is no transmon) but kept so the notebook's import line
# resolves.
L_PROBE_NH = L_J_NH
CJ_PROBE_FF = CJ_J_FF
L_TRANSMON_NH = L_J_NH
CJ_TRANSMON_FF = CJ_J_FF
QUBIT_COMPONENT_NAME = ARRAY_COMPONENT_NAME

# ---- chip geometry ----
# size_x/size_y/center_x/center_y are auto-computed from the drawn geometry's
# bounding box (see _design_bounds_mm / CHIP_MARGIN_FACTOR) unless explicitly
# overridden. size_z is fixed. The notebook itself leaves the chip at the
# DesignPlanar default; here it is sized to the geometry so the mesh domain is
# not needlessly large.
CHIP_SIZE_Z_UM = 500
# chip footprint = design bounding box * this factor. 2.5 (was 4.0): the
# geometry is ~1.3mm and a 4x margin gave a 4.7x3.7mm chip whose far vacuum/
# substrate volume dominated the FEM mesh and spawned box modes near the
# resonator. 2.5 keeps every chip edge >500um (>>25 trace widths) from any
# trace - no effect on the resonator frequency - while cutting ~40% of the
# mesh volume and pushing box modes well above the resonator.
CHIP_MARGIN_FACTOR = 2.5


def identify_resonator_mode(mat_mode_port, freqs_hz, phonon_abs_sum_thresh=1e-2,
                            min_ghz=None):
    """Pick the readout-resonator eigenmode out of a solve that also contains
    SQUID-array phonon modes.

    The discriminator is the sum of ABSOLUTE per-junction participations,
    ``sum_j |p_{m,j}|``, over the array's inductive ports:

      * a phonon stores ~all of its inductive energy in the array, so
        ``sum|p_j| ~= 1`` (individual p_j alternate sign, so the *signed*
        sum can still be tiny - which is why signed ``p_array`` fails here);
      * a mode on the CPW (the readout resonator) or a package/box mode
        barely touches the array, so ``sum|p_j| << 1e-3``.

    Among the non-array modes the readout resonator is the lowest in
    frequency (its lambda/4 fundamental; box modes and CPW harmonics are
    higher). Returns ``(index, info)``; raises if no non-array mode is found.
    """
    abs_sum = np.abs(mat_mode_port).sum(axis=1)
    cand = [i for i in range(len(freqs_hz))
            if abs_sum[i] < phonon_abs_sum_thresh
            and (min_ghz is None or freqs_hz[i] / 1e9 >= min_ghz)]
    if not cand:
        raise RuntimeError(
            f"identify_resonator_mode: no non-array mode (min sum|p_j| = "
            f"{abs_sum.min():.2e} >= threshold {phonon_abs_sum_thresh}). "
            f"Widen the solve window or raise the threshold.")
    idx = min(cand, key=lambda i: freqs_hz[i])
    return idx, {
        "abs_sum_p": float(abs_sum[idx]),
        "phonon_abs_sum_thresh": phonon_abs_sum_thresh,
        "n_nonarray_candidates": len(cand),
        "nonarray_candidate_ghz": sorted(float(freqs_hz[i] / 1e9) for i in cand),
        "abs_sum_p_all_modes": [float(x) for x in abs_sum],
    }


def _design_bounds_mm(design):
    """Bounding box (minx, miny, maxx, maxy), in mm, across every drawn
    qgeometry table. Must be called AFTER all components are built."""
    all_bounds = []
    for table_name, table in design.qgeometry.tables.items():
        if len(table) == 0:
            continue
        all_bounds.append(table.total_bounds)  # [minx, miny, maxx, maxy], mm
    all_bounds = np.array(all_bounds)
    return (all_bounds[:, 0].min(), all_bounds[:, 1].min(),
            all_bounds[:, 2].max(), all_bounds[:, 3].max())


def build_design(n_squids=N_SQUIDS, L_j_nH=L_J_NH, res_length_mm=RES_LENGTH_MM,
                 squid_gap_um=SQUID_GAP_UM, arm_width_um=ARM_WIDTH_UM,
                 josephson_gap_um=JOSEPHSON_GAP_UM, segment_length_um=SEGMENT_LENGTH_UM,
                 arm_semi_length_um=ARM_SEMI_LENGTH_UM, line_width_um=LINE_WIDTH_UM,
                 ground_gap_um=GROUND_GAP_UM, jj_width_um=JJ_WIDTH_UM, cj_j_fF=CJ_J_FF,
                 array_pos_x_um=ARRAY_POS_X_UM, array_pos_y_um=ARRAY_POS_Y_UM,
                 chip_size_x_mm=None, chip_size_y_mm=None,
                 chip_center_x_mm=None, chip_center_y_mm=None):
    """
    Build the SQUID-array + meandered readout-resonator design drawn in
    build_design.ipynb.

    The array is a series chain of ``n_squids`` SquidArray loops; every one of
    its ``2 * n_squids`` junctions carries inductance ``L_j_nH``. The resonator
    runs between an OpenToGround and a ShortToGround termination and is not
    galvanically connected to the array - coupling is purely the resonator
    current threading flux through the loops.

    The chip is sized AFTER the geometry is drawn, from the geometry's own
    bounding box * CHIP_MARGIN_FACTOR, centered on that bounding box. Any of the
    chip_* kwargs overrides the corresponding auto-computed value.
    """
    design = designs.DesignPlanar({}, overwrite_enabled=True)

    design.variables['cpw_width'] = '20um'
    design.variables['cpw_gap'] = '12.25um'

    # --- SQUID array (series chain = tunable superinductor) ---

    SquidArray(design, ARRAY_COMPONENT_NAME, options={
        'pos_x': f'{array_pos_x_um}um',
        'pos_y': f'{array_pos_y_um}um',
        'orientation': f'{ARRAY_ORIENTATION_DEG}',
        'n_squids': f'{int(n_squids)}',
        'plate1_width': f'{arm_width_um}um',
        'squid_gap': f'{squid_gap_um}um',
        'segment_a_length': f'{arm_semi_length_um}um',
        'segment_a_width': f'{arm_width_um}um',
        'JJ_gap': f'{josephson_gap_um}um',
        'segment_b_length': f'{arm_semi_length_um}um',
        'segment_b_width': f'{arm_width_um}um',
        'segment_c_width': f'{arm_width_um}um',
        'segment_d_length': f'{segment_length_um}um',
        'segment_d_width': f'{line_width_um}um',
        'jj_width': f'{jj_width_um}um',
        'L_j': f'{L_j_nH}nH',
        'C_j': f'{cj_j_fF}fF',
        'ground_gap': f'{ground_gap_um}um',
    })

    # --- Readout resonator terminations ---

    OpenToGround(design, 'launch_point1', options=dict(
        pos_x='-275um', pos_y='-850um', orientation='0',
        termination_gap='cpw_gap', gap='cpw_gap', width='cpw_width'))

    ShortToGround(design, 'launch_point2', options=dict(
        pos_x='-150um', pos_y='50um', orientation='0',
        termination_gap='cpw_gap', gap='cpw_gap', width='cpw_width'))

    jogs = OrderedDict()
    jogs[0] = ["L", "100um"]
    jogs[1] = ["L", "700um"]

    RouteMeander(design, 'readout_1', options=dict(
        pin_inputs=dict(
            start_pin=dict(component='launch_point1', pin='open'),
            end_pin=dict(component='launch_point2', pin='short'),
        ),
        fillet='40um',
        lead=dict(
            start_straight='10um',
            end_straight='1180um',
            end_jogged_extension=jogs,
        ),
        total_length=f'{res_length_mm}mm',
    ))

    # --- Size the chip from the drawn geometry's own bounding box ---

    minx, miny, maxx, maxy = _design_bounds_mm(design)
    bounds_center_x, bounds_center_y = (minx + maxx) / 2, (miny + maxy) / 2
    bounds_size_x, bounds_size_y = (maxx - minx), (maxy - miny)

    chip_size_x_mm = bounds_size_x * CHIP_MARGIN_FACTOR if chip_size_x_mm is None else chip_size_x_mm
    chip_size_y_mm = bounds_size_y * CHIP_MARGIN_FACTOR if chip_size_y_mm is None else chip_size_y_mm
    chip_center_x_mm = bounds_center_x if chip_center_x_mm is None else chip_center_x_mm
    chip_center_y_mm = bounds_center_y if chip_center_y_mm is None else chip_center_y_mm

    design.chips.main.size.size_x = f'{chip_size_x_mm}mm'
    design.chips.main.size.size_y = f'{chip_size_y_mm}mm'
    design.chips.main.size.size_z = f'{CHIP_SIZE_Z_UM}um'
    design.chips.main.size.center_x = f'{chip_center_x_mm}mm'
    design.chips.main.size.center_y = f'{chip_center_y_mm}mm'

    return design


if __name__ == "__main__":
    print(f"N_SQUIDS={N_SQUIDS}, L_J_NH={L_J_NH}, RES_LENGTH_MM={RES_LENGTH_MM}")
    d = build_design()
    d.rebuild()
    n_jj = len(d.qgeometry.tables['junction'].loc[
        d.qgeometry.tables['junction']['component'] == d.components[ARRAY_COMPONENT_NAME].id])
    mnx, mny, mxx, mxy = _design_bounds_mm(d)
    print("Design built OK. Components:", list(d.components.keys()))
    print(f"Array junction rows: {n_jj} (expect {JJ_PER_LOOP * N_SQUIDS})")
    print(f"Design bounds (mm): X=[{mnx:.4f},{mxx:.4f}]  Y=[{mny:.4f},{mxy:.4f}]")
    print(f"Chip size (mm): {d.chips.main.size.size_x} x {d.chips.main.size.size_y}, "
          f"centered at ({d.chips.main.size.center_x}, {d.chips.main.size.center_y})")
