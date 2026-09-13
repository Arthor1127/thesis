from collections import OrderedDict

import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.qubits.squid_array import SquidArray
from qiskit_metal.qlibrary.sample_shapes.rect_cutout import RectCutout
from qiskit_metal.qlibrary.sample_shapes.rectangle import Rectangle

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
GROUND_GAP_UM = 3.0        # -> ground_gap       (SquidArray's own conformal JJ_gap clearance)
ARRAY_GROUND_MARGIN_UM = 6.0  # clean rectangular ground pocket: array envelope + this margin
                             # (overrides SquidArray's buffer(ground_gap) pocket, which pinches
                             #  in over the thin inter-loop leads and floods ground between loops).
                             # Kept < route_clearance - 22.25 um (CPW half-gap) so a ground strip
                             # still survives between the hairpin and the array at the tightest
                             # route_clearance in the sweep (30 um).

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

# ---- readout line ----
# Vestigial: the readout line is now a fixed-geometry hairpin (RouteAnchors), not a
# length-targeted meander. Kept only so build_design.ipynb's `from build_design
# import ... RES_LENGTH_MM` still resolves; `res_length_mm` is ignored.
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
# For the magnetostatic inductance sim, pass chip_margin_um instead: the hairpin
# field is a localised dipole (opposite currents in the two passes) that decays
# as 1/r^2, so a few hundred um of ground plane past the geometry is plenty and
# an absolute margin beats a 2.5x multiplier on the ~1.2mm bbox.
CHIP_MARGIN_FACTOR = 2.5

# ---- readout line: hairpin around the array ----
# gap from the readout-line centreline to the nearest array metal, held equal on
# the top / bottom / left-elbow / termination sides
ROUTE_CLEARANCE_UM = 35.0
ROUTE_FILLET_UM = 20.0    # keep < ROUTE_CLEARANCE_UM so the elbow arc clears the array


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


def loop_aperture_polys(n_squids=N_SQUIDS, squid_gap_um=SQUID_GAP_UM,
                        arm_width_um=ARM_WIDTH_UM, arm_semi_length_um=ARM_SEMI_LENGTH_UM,
                        josephson_gap_um=JOSEPHSON_GAP_UM, segment_length_um=SEGMENT_LENGTH_UM,
                        array_pos_x_um=ARRAY_POS_X_UM, array_pos_y_um=ARRAY_POS_Y_UM,
                        grow_um=0.3):
    """One shapely box per SQUID-loop opening, **in metres**, computed from the
    array parameters (no need to build the SquidArray component). Geometry matches
    ``squid_array.SquidArray._loop``: aperture i is ``loop_w`` wide x ``squid_gap``
    tall, its left edge at plate1's right edge.

    ``grow_um`` (default 0.3) expands every box outward on all sides so its edges
    do NOT coincide with the array-metal edges -- an exact coincidence makes GMSH
    build degenerate zero-length edges when it partitions the dielectric gap by
    the integration surface ("Could not create line"). The extra sliver sits over
    PEC metal where B_z ~ 0, so it adds ~no flux and, being a fixed 0.3 um, adds
    the same ~nothing regardless of aperture size (no size-dependent bias).

    Returns a ``shapely.MultiPolygon``. Use it as the integration area and to size
    the ground-plane window when the array metal is left out of the EM model.
    """
    import shapely
    full_w = (arm_width_um + 2 * arm_semi_length_um + josephson_gap_um
              + arm_width_um + segment_length_um)
    loop_w = 2 * arm_semi_length_um + josephson_gap_um
    g = grow_um
    rects = []
    for i in range(int(n_squids)):
        x_left = array_pos_x_um + i * full_w + 0.5 * arm_width_um
        rects.append(shapely.geometry.box(
            (x_left - g) * 1e-6, (array_pos_y_um - 0.5 * squid_gap_um - g) * 1e-6,
            (x_left + loop_w + g) * 1e-6, (array_pos_y_um + 0.5 * squid_gap_um + g) * 1e-6))
    return shapely.MultiPolygon(rects)


def geometry_anchors(n_squids=N_SQUIDS, squid_gap_um=SQUID_GAP_UM,
                     arm_width_um=ARM_WIDTH_UM, arm_semi_length_um=ARM_SEMI_LENGTH_UM,
                     josephson_gap_um=JOSEPHSON_GAP_UM, segment_length_um=SEGMENT_LENGTH_UM,
                     array_pos_x_um=ARRAY_POS_X_UM, array_pos_y_um=ARRAY_POS_Y_UM,
                     route_clearance_um=ROUTE_CLEARANCE_UM, term_x_um=None):
    """The array/hairpin layout anchors ``build_design`` computes, as a dict (um)
    -- single source of truth, so post-processing scripts (e.g. eigenmode field
    extraction) don't re-derive these formulas by hand. Matches ``build_design``
    exactly; keep the two in sync if either changes.

    Keys: ``full_w`` (loop pitch), ``arr_left``/``arr_right``/``arr_half_h``
    (array metal bounding box), ``loop_x0`` (array of each loop's aperture left
    edge x, length ``n_squids``), ``pass_y_top``/``pass_y_bot`` (hairpin pass
    centrelines), ``elbow_x`` (hairpin U-turn x), ``term_x`` (both terminations'
    x, i.e. the OpenToGround/ShortToGround pin location).
    """
    full_w = (arm_width_um + 2 * arm_semi_length_um + josephson_gap_um
             + arm_width_um + segment_length_um)
    arr_left = array_pos_x_um - 0.5 * arm_width_um
    arr_right = array_pos_x_um + n_squids * full_w - 0.5 * arm_width_um
    arr_half_h = 0.5 * squid_gap_um + arm_width_um
    loop_x0 = array_pos_x_um + np.arange(int(n_squids)) * full_w + 0.5 * arm_width_um
    pass_y_top = array_pos_y_um + arr_half_h + route_clearance_um
    pass_y_bot = array_pos_y_um - arr_half_h - route_clearance_um
    elbow_x = arr_left - route_clearance_um
    if term_x_um is None:
        term_x_um = arr_right + route_clearance_um
    return dict(full_w=full_w, arr_left=arr_left, arr_right=arr_right,
               arr_half_h=arr_half_h, loop_x0=loop_x0, pass_y_top=pass_y_top,
               pass_y_bot=pass_y_bot, elbow_x=elbow_x, term_x=term_x_um)


def build_design(n_squids=N_SQUIDS, L_j_nH=L_J_NH, res_length_mm=RES_LENGTH_MM,
                 squid_gap_um=SQUID_GAP_UM, arm_width_um=ARM_WIDTH_UM,
                 josephson_gap_um=JOSEPHSON_GAP_UM, segment_length_um=SEGMENT_LENGTH_UM,
                 arm_semi_length_um=ARM_SEMI_LENGTH_UM, line_width_um=LINE_WIDTH_UM,
                 ground_gap_um=GROUND_GAP_UM, array_ground_margin_um=ARRAY_GROUND_MARGIN_UM,
                 jj_width_um=JJ_WIDTH_UM, cj_j_fF=CJ_J_FF,
                 array_pos_x_um=ARRAY_POS_X_UM, array_pos_y_um=ARRAY_POS_Y_UM,
                 route_clearance_um=ROUTE_CLEARANCE_UM, route_fillet_um=ROUTE_FILLET_UM,
                 term_x_um=None, cpw_width_um=20.0, cpw_gap_um=12.25,
                 chip_size_x_mm=None, chip_size_y_mm=None,
                 chip_center_x_mm=None, chip_center_y_mm=None,
                 chip_margin_um=None, chip_size_z_um=CHIP_SIZE_Z_UM,
                 include_array_metal=True, ground_near_end=False):
    """
    Build the SQUID-array + readout-line design drawn in build_design.ipynb.

    The array is a series chain of ``n_squids`` SquidArray loops; every one of
    its ``2 * n_squids`` junctions carries inductance ``L_j_nH``. The readout line
    is a hairpin (RouteAnchors, no prescribed length) that runs past the top of
    the array, U-turns at the left end, and comes back past the bottom, between an
    OpenToGround and a ShortToGround termination. It is not galvanically connected
    to the array - coupling is purely the line current threading flux through the
    loops.

    ``route_clearance_um`` is the gap from the readout-line centreline to the
    nearest array metal, applied equally on the top, bottom, left-elbow and
    termination sides (``term_x_um`` overrides the last). ``route_fillet_um`` must
    stay below ``route_clearance_um`` so the elbow arc clears the array.

    The chip is sized AFTER the geometry is drawn. By default it is the geometry
    bounding box * ``CHIP_MARGIN_FACTOR``; pass ``chip_margin_um`` to use an
    absolute margin instead (bbox + 2 * margin), which is what the inductance sim
    wants to keep the FEM domain small. Any ``chip_*`` kwarg overrides the
    corresponding auto-computed value; ``chip_size_z_um`` sets the substrate
    thickness (magnetostatically inert - silicon has mu_r = 1 - so it only trades
    mesh volume).

    ``include_array_metal`` (default True) draws the real ``SquidArray``. Set it
    False for the inductance sim to leave the array metal *out* of the EM model:
    the ~23 open SQUID loops are topological handles in the field domain that
    wreck the AMS solver's convergence (~500 iterations, forces GMRES). Without
    the metal the loop openings are just vacuum apertures, a plain
    ``RectCutout`` opens the ground plane over them, and the solve converges in
    ~100 CG iterations. Use ``loop_aperture_polys(...)`` (same signature args) for
    the integration area. Whether the two agree is a question the sim answers.

    ``ground_near_end`` (default False) adds a positive-metal bridge tying the
    array's left end (pin 'a', nearest the resonator elbow) to the ground plane,
    for the one-off "does grounding the near end change the flux" comparison.

    ``array_ground_margin_um`` (default ARRAY_GROUND_MARGIN_UM) cuts a clean
    rectangular ground pocket over the array envelope + this margin. SquidArray's
    own pocket is ``all_metal.buffer(ground_gap)``, which follows the metal
    outline and so pinches down to ~``line_width`` over each inter-loop lead,
    letting the ground plane flood into the notches between loops -- grounded PEC
    tabs a few um from every loop, whose size tracks ``segment_length`` and which
    focus/screen the hairpin field right at the loop apertures. The rectangular
    pocket removes them. Set to 0/None to keep SquidArray's conformal pocket.
    """
    design = designs.DesignPlanar({}, overwrite_enabled=True)

    design.variables['cpw_width'] = f'{cpw_width_um}um'
    design.variables['cpw_gap'] = f'{cpw_gap_um}um'

    # one loop pitch / aperture, matching squid_array.SquidArray.make()
    full_w_um = (arm_width_um + 2 * arm_semi_length_um + josephson_gap_um
                 + arm_width_um + segment_length_um)
    arr_left_um = array_pos_x_um - 0.5 * arm_width_um
    arr_right_um = array_pos_x_um + n_squids * full_w_um - 0.5 * arm_width_um
    arr_half_h_um = 0.5 * squid_gap_um + arm_width_um     # loop centre -> top metal edge

    if include_array_metal:
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
        # Clean rectangular ground pocket over the whole array envelope. Unions
        # with (and dominates) SquidArray's conformal buffer(ground_gap) pocket,
        # which otherwise pinches down to ~line_width over each inter-loop lead
        # and lets the ground plane flood into the notches between loops (a
        # segment_length-dependent screen right in the hairpin->loop field path,
        # and ground_gap-scale re-entrant corners near the fab resolution limit).
        # A straight-edged pocket with a uniform margin fixes both.
        if array_ground_margin_um and array_ground_margin_um > 0:
            RectCutout(design, ARRAY_COMPONENT_NAME + '_pocket', options=dict(
                pos_x=f'{0.5 * (arr_left_um + arr_right_um)}um',
                pos_y=f'{array_pos_y_um}um',
                width=f'{(arr_right_um - arr_left_um) + 2 * array_ground_margin_um}um',
                height=f'{2 * arr_half_h_um + 2 * array_ground_margin_um}um'))

        if ground_near_end:
            # tie the array's LEFT end (pin 'a', the one nearest the resonator
            # elbow) to the ground plane: a positive-metal bridge spanning the
            # pocket. Positive polys survive the pocket subtraction. Reach 3 um
            # past the (rectangular or conformal) pocket edge into solid ground.
            pocket_edge_x = arr_left_um - (array_ground_margin_um
                                          if array_ground_margin_um and array_ground_margin_um > 0
                                          else ground_gap_um)
            bridge_x0 = pocket_edge_x - 3.0                   # 3 um into ground
            bridge_x1 = arr_left_um + 2.0                     # 2 um into plate1
            Rectangle(design, ARRAY_COMPONENT_NAME + '_gnd_bridge', options=dict(
                pos_x=f'{0.5 * (bridge_x0 + bridge_x1)}um',
                pos_y=f'{array_pos_y_um}um',
                width=f'{bridge_x1 - bridge_x0}um',
                height=f'{2 * arm_width_um}um'))
    else:
        # array metal omitted: open the ground plane over the loop band instead
        # (matches the array's own ground_pocket = metal bbox + ground_gap)
        RectCutout(design, ARRAY_COMPONENT_NAME + '_window', options=dict(
            pos_x=f'{0.5 * (arr_left_um + arr_right_um)}um',
            pos_y=f'{array_pos_y_um}um',
            width=f'{(arr_right_um - arr_left_um) + 2 * ground_gap_um}um',
            height=f'{2 * arr_half_h_um + 2 * ground_gap_um}um'))

    # --- Readout line: a hairpin straddling the array (length just follows) ---
    #
    # RouteMeander needs total_length and adds wiggles to hit it; here the wire is
    # a plain U whose length falls out of the geometry. It runs past the top of the
    # array, U-turns at the left end, and comes back past the bottom. Every side is
    # route_clearance_um from the nearest array metal. Current enters at
    # launch_point1 (open) and returns to ground at launch_point2 (short); its
    # field threads every SQUID loop. (arr_* and full_w_um computed above.)

    pass_y_top_um = array_pos_y_um + arr_half_h_um + route_clearance_um
    pass_y_bot_um = array_pos_y_um - arr_half_h_um - route_clearance_um
    elbow_x_um = arr_left_um - route_clearance_um
    if term_x_um is None:
        term_x_um = arr_right_um + route_clearance_um

    OpenToGround(design, 'launch_point1', options=dict(
        pos_x=f'{term_x_um}um', pos_y=f'{pass_y_top_um}um', orientation='0',
        termination_gap='cpw_gap', gap='cpw_gap', width='cpw_width'))

    ShortToGround(design, 'launch_point2', options=dict(
        pos_x=f'{term_x_um}um', pos_y=f'{pass_y_bot_um}um', orientation='0',
        termination_gap='cpw_gap', gap='cpw_gap', width='cpw_width'))

    anchors = OrderedDict()
    anchors[0] = np.array([elbow_x_um, pass_y_top_um]) * 1e-3    # mm
    anchors[1] = np.array([elbow_x_um, pass_y_bot_um]) * 1e-3

    RouteAnchors(design, 'readout_1', options=dict(
        pin_inputs=dict(
            start_pin=dict(component='launch_point1', pin='open'),
            end_pin=dict(component='launch_point2', pin='short'),
        ),
        fillet=f'{route_fillet_um}um',
        trace_width='cpw_width',
        trace_gap='cpw_gap',
        lead=dict(start_straight='30um', end_straight='30um'),
        anchors=anchors,
    ))

    # --- Size the chip from the drawn geometry's own bounding box ---

    minx, miny, maxx, maxy = _design_bounds_mm(design)
    bounds_center_x, bounds_center_y = (minx + maxx) / 2, (miny + maxy) / 2
    bounds_size_x, bounds_size_y = (maxx - minx), (maxy - miny)

    if chip_margin_um is not None:
        m = chip_margin_um * 1e-3
        chip_size_x_mm = (bounds_size_x + 2 * m) if chip_size_x_mm is None else chip_size_x_mm
        chip_size_y_mm = (bounds_size_y + 2 * m) if chip_size_y_mm is None else chip_size_y_mm
    else:
        chip_size_x_mm = bounds_size_x * CHIP_MARGIN_FACTOR if chip_size_x_mm is None else chip_size_x_mm
        chip_size_y_mm = bounds_size_y * CHIP_MARGIN_FACTOR if chip_size_y_mm is None else chip_size_y_mm
    chip_center_x_mm = bounds_center_x if chip_center_x_mm is None else chip_center_x_mm
    chip_center_y_mm = bounds_center_y if chip_center_y_mm is None else chip_center_y_mm

    design.chips.main.size.size_x = f'{chip_size_x_mm}mm'
    design.chips.main.size.size_y = f'{chip_size_y_mm}mm'
    design.chips.main.size.size_z = f'{chip_size_z_um}um'
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
