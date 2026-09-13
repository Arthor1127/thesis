"""Bare Josephson-junction array for an eigenmode ("normal modes") solve.

This is the geometry of ``josephson_array_blake.ipynb`` (a series chain of
identical Dolan-bridge junctions) reduced to the simulation-oriented
``JosephsonArray`` qiskit-metal component: ``n_junctions`` rectangular islands
separated by ``pad_gap`` vacuum breaks, one lumped-inductor ``junction`` line per
break.

There is NO readout resonator here and NO EPR / participation analysis - the
only thing wanted is the set of normal-mode frequencies of the array itself
(its standing-wave plasma modes), so the chain sits alone in a large
ground-plane pocket and every junction gap is bridged by a lumped ``L_j``/``C_j``
port in ``eigenmodes.py``.

``JosephsonArray`` lives in the ``qcg-quantum-design`` env at
``qiskit_metal/qlibrary/qubits/josephson_array.py`` (added alongside
``squid_array.py``); it must be present in the cluster env too.
"""

import numpy as np
from qiskit_metal import designs
from qiskit_metal.qlibrary.qubits.josephson_array import JosephsonArray

# ---- array: junction electrical parameters ----
N_JUNCTIONS = 40
L_J_NH = 1.2         # inductance of EVERY junction (JosephsonArray.L_j)
CJ_J_FF = 2.0        # capacitance of every junction (JosephsonArray.C_j)
JJ_WIDTH_UM = 0.5    # transverse width of each junction bridge

# ---- array: island geometry (um) ----
PAD_LENGTH_UM = 8.0    # island length along the chain (-> pad_length)
PAD_WIDTH_UM = 6.0     # island transverse width          (-> pad_width)
PAD_GAP_UM = 0.15      # vacuum break each junction bridges (-> pad_gap)

# ---- ground-plane pocket ----
# "bigger cutout of the ground plane": a 50um trench around the whole chain,
# well clear of the ~6-8um island features, so the array's plasma modes are
# not loaded / shorted by the ground plane. (SquidArray runs used only 3um.)
GROUND_GAP_UM = 50.0

# ---- placement ----
ARRAY_POS_X_UM = 0.0
ARRAY_POS_Y_UM = 0.0
ARRAY_ORIENTATION_DEG = 0.0
ARRAY_COMPONENT_NAME = 'array_1'

# ---- chip / air-box ----
# The array is long and thin, so the chip footprint is the geometry bounding
# box plus a fixed margin on every side (an additive margin, not a x-factor:
# a multiplicative factor would leave the thin y direction too tight and blow
# the long x direction up). size_z is fixed.
CHIP_MARGIN_UM = 300.0
CHIP_SIZE_Z_UM = 500.0


def _design_bounds_mm(design):
    """Bounding box (minx, miny, maxx, maxy) in mm across every qgeometry
    table. Call AFTER all components are built."""
    all_bounds = []
    for _, table in design.qgeometry.tables.items():
        if len(table) == 0:
            continue
        all_bounds.append(table.total_bounds)  # [minx, miny, maxx, maxy], mm
    all_bounds = np.array(all_bounds)
    return (all_bounds[:, 0].min(), all_bounds[:, 1].min(),
            all_bounds[:, 2].max(), all_bounds[:, 3].max())


def build_design(n_junctions=N_JUNCTIONS, L_j_nH=L_J_NH, cj_j_fF=CJ_J_FF,
                 jj_width_um=JJ_WIDTH_UM, pad_length_um=PAD_LENGTH_UM,
                 pad_width_um=PAD_WIDTH_UM, pad_gap_um=PAD_GAP_UM,
                 ground_gap_um=GROUND_GAP_UM,
                 array_pos_x_um=ARRAY_POS_X_UM, array_pos_y_um=ARRAY_POS_Y_UM,
                 chip_margin_um=CHIP_MARGIN_UM,
                 chip_size_x_mm=None, chip_size_y_mm=None,
                 chip_center_x_mm=None, chip_center_y_mm=None):
    """Build the bare JJ-array design.

    ``n_junctions`` junctions -> ``n_junctions`` lumped inductive ports (created
    in ``eigenmodes.py``). The chip is sized from the drawn geometry's bounding
    box + ``chip_margin_um`` on each side unless a ``chip_*`` kwarg overrides it.
    """
    design = designs.DesignPlanar({}, overwrite_enabled=True)

    JosephsonArray(design, ARRAY_COMPONENT_NAME, options={
        'pos_x': f'{array_pos_x_um}um',
        'pos_y': f'{array_pos_y_um}um',
        'orientation': f'{ARRAY_ORIENTATION_DEG}',
        'n_junctions': f'{int(n_junctions)}',
        'pad_length': f'{pad_length_um}um',
        'pad_width': f'{pad_width_um}um',
        'pad_gap': f'{pad_gap_um}um',
        'jj_width': f'{jj_width_um}um',
        'L_j': f'{L_j_nH}nH',
        'C_j': f'{cj_j_fF}fF',
        'ground_gap': f'{ground_gap_um}um',
    })

    design.rebuild()

    # --- size the chip from the drawn geometry's own bounding box ---
    minx, miny, maxx, maxy = _design_bounds_mm(design)
    margin_mm = chip_margin_um * 1e-3
    bounds_center_x, bounds_center_y = (minx + maxx) / 2, (miny + maxy) / 2

    chip_size_x_mm = (maxx - minx) + 2 * margin_mm if chip_size_x_mm is None else chip_size_x_mm
    chip_size_y_mm = (maxy - miny) + 2 * margin_mm if chip_size_y_mm is None else chip_size_y_mm
    chip_center_x_mm = bounds_center_x if chip_center_x_mm is None else chip_center_x_mm
    chip_center_y_mm = bounds_center_y if chip_center_y_mm is None else chip_center_y_mm

    design.chips.main.size.size_x = f'{chip_size_x_mm}mm'
    design.chips.main.size.size_y = f'{chip_size_y_mm}mm'
    design.chips.main.size.size_z = f'{CHIP_SIZE_Z_UM}um'
    design.chips.main.size.center_x = f'{chip_center_x_mm}mm'
    design.chips.main.size.center_y = f'{chip_center_y_mm}mm'

    return design


if __name__ == "__main__":
    for n in (10, 23, 40, 100):
        d = build_design(n_junctions=n)
        n_jj = len(d.qgeometry.tables['junction'])
        mnx, mny, mxx, mxy = _design_bounds_mm(d)
        print(f"n_junctions={n:3d}: junction rows={n_jj:3d}  "
              f"bbox(mm) X=[{mnx:.4f},{mxx:.4f}] Y=[{mny:.4f},{mxy:.4f}]  "
              f"chip={d.chips.main.size.size_x} x {d.chips.main.size.size_y}")
