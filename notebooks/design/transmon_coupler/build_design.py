from collections import OrderedDict

import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.qubits.circle_transmon import CircTransmon
from qiskit_metal.qlibrary.couplers.cap_n_interdigital_round import CapNInterdigitalRound

# ---------------------------------------------------------------------------
# Two grounded circular transmons (CircTransmon) sharing a single fixed
# capacitive coupler. The coupler (CapNInterdigitalRound) is an interdigitated
# finger capacitor with fully rounded metal - no sharp corners - built so that
# each of its two combs is a direct galvanic extension of one transmon island.
# The coupler is told the pad geometry (island_radius / island_gap /
# island_center_offset), so each comb's lead flares out and blends tangentially
# into the round pad, and its ground moat merges into each transmon pocket with
# the cusps filleted. The two combs never touch, so Q1 and Q2 stay on separate
# nets and the inter-comb capacitance is the coupling capacitance C_c.
#
# No CPW / claws / routing. The chip is sized AFTER the geometry is drawn, from
# the geometry's own bounding box.
# ---------------------------------------------------------------------------

# ---- transmon junctions (one JJ to ground per qubit) ----
L_Q1_NH = 10.0
L_Q2_NH = 10.0
CJ_Q1_FF = 2.0
CJ_Q2_FF = 2.0
JJ_WIDTH_UM = 5.0
JJ_SIM_GAP_UM = 3.0

# ---- transmon pad geometry (um) ----
PAD_RADIUS_UM = 136.0
PAD_GAP_UM = 20.0
JJ_ANGLE_DEG = 270            # JJ points south, away from the east/west coupler

# ---- coupler: interdigital capacitor geometry (um) ----
CPL_FINGER_COUNT = 8
CPL_FINGER_WIDTH_UM = 8.0
CPL_FINGER_LENGTH_UM = 80.0
CPL_FINGER_GAP_UM = 4.0
CPL_FINGER_TIP_GAP_UM = 6.0
CPL_SPINE_WIDTH_UM = 12.0
CPL_FILLET_RADIUS_UM = 6.0
CPL_CAP_GAP_GROUND_UM = 3.0
CPL_ETCH_CORNER_RADIUS_UM = 6.0     # corner radius of the rectangular ground moat
CPL_LEAD_LENGTH_UM = 40.0            # spine-to-pad gap the neck + flare span
CPL_ROUND_RESOLUTION = 32

# circular-island blending: the coupler is told the transmon pad geometry so it
# leaves each pad as a narrow neck that flares tangentially into the round pad
# (large lead_blend_radius = sweeping trumpet), and so its moat merges into each
# transmon pocket with the cusps filleted (no sharp points at the junction).
CPL_LEAD_WIDTH_UM = 12.0             # neck width where it crosses the pad edge
CPL_LEAD_FLARE_WIDTH_UM = 44.0       # neck width where it meets the comb spine
CPL_LEAD_BLEND_RADIUS_UM = 40.0      # tangent-flare radius at the pad
CPL_ETCH_BLEND_RADIUS_UM = 10.0

# how far each transmon pad edge overlaps the coupler lead-stub tip
ABUT_OVERLAP_UM = 4.0

# ---- component names ----
Q1_NAME = 'Q1'
Q2_NAME = 'Q2'
COUPLER_NAME = 'coupler'

# SQDMetal's create_port_JosephsonJunction addresses junctions positionally
# (insertion order into design.qgeometry.tables['junction']). Each CircTransmon
# emits exactly one junction, so Q1 -> junction index 0, Q2 -> junction index 1
# (Q1 is built first below). Per component the JJ is index 0.
JJ_INDEX = 0

# ---- chip geometry ----
CHIP_SIZE_Z_UM = 500
CHIP_MARGIN_FACTOR = 3.0       # chip footprint = design bounding box * this factor


def _coupler_half_span_um():
    """Distance (um) from the coupler centre to a lead-stub tip, in the local
    frame: half the finger channel + one spine + the lead stub."""
    half_h = 0.5 * (CPL_FINGER_LENGTH_UM + CPL_FINGER_TIP_GAP_UM)
    return half_h + CPL_SPINE_WIDTH_UM + CPL_LEAD_LENGTH_UM


def _transmon_offset_um():
    """|pos_x| of each transmon so its pad edge overlaps the coupler stub tip by
    ABUT_OVERLAP_UM."""
    return PAD_RADIUS_UM + _coupler_half_span_um() - ABUT_OVERLAP_UM


def _design_bounds_mm(design):
    """Bounding box (minx, miny, maxx, maxy), in mm, across every drawn
    qgeometry table. Must be called AFTER all components are built."""
    all_bounds = []
    for _, table in design.qgeometry.tables.items():
        if len(table) == 0:
            continue
        all_bounds.append(table.total_bounds)  # [minx, miny, maxx, maxy], mm
    all_bounds = np.array(all_bounds)
    return (all_bounds[:, 0].min(), all_bounds[:, 1].min(),
            all_bounds[:, 2].max(), all_bounds[:, 3].max())


def _transmon_options(pos_x_um, L_nH, Cj_fF):
    return Dict(
        pos_x=f'{pos_x_um}um', pos_y='0um', orientation='0', chip='main',
        cpw_width='cpw_width', cpw_gap='cpw_gap',
        pad_radius=f'{PAD_RADIUS_UM}um',
        pad_gap=f'{PAD_GAP_UM}um',
        jj_options=Dict(
            jj_width=f'{JJ_WIDTH_UM}um',
            jj_angle=f'{JJ_ANGLE_DEG}',
            L_j=f'{L_nH}nH',
            C_j=f'{Cj_fF}fF',
            export_mask=False,
            jj_sim_gap=f'{JJ_SIM_GAP_UM}um',
        ),
    )


def _coupler_options():
    return Dict(
        pos_x='0um', pos_y='0um', orientation='90', chip='main',
        finger_count=f'{CPL_FINGER_COUNT}',
        finger_width=f'{CPL_FINGER_WIDTH_UM}um',
        finger_length=f'{CPL_FINGER_LENGTH_UM}um',
        finger_gap=f'{CPL_FINGER_GAP_UM}um',
        finger_tip_gap=f'{CPL_FINGER_TIP_GAP_UM}um',
        spine_width=f'{CPL_SPINE_WIDTH_UM}um',
        fillet_radius=f'{CPL_FILLET_RADIUS_UM}um',
        cap_gap_ground=f'{CPL_CAP_GAP_GROUND_UM}um',
        etch_corner_radius=f'{CPL_ETCH_CORNER_RADIUS_UM}um',
        lead_length=f'{CPL_LEAD_LENGTH_UM}um',
        lead_width=f'{CPL_LEAD_WIDTH_UM}um',
        round_resolution=f'{CPL_ROUND_RESOLUTION}',
        island_radius=f'{PAD_RADIUS_UM}um',
        island_gap=f'{PAD_GAP_UM}um',
        island_center_offset=f'{_transmon_offset_um()}um',
        lead_flare_width=f'{CPL_LEAD_FLARE_WIDTH_UM}um',
        lead_blend_radius=f'{CPL_LEAD_BLEND_RADIUS_UM}um',
        etch_blend_radius=f'{CPL_ETCH_BLEND_RADIUS_UM}um',
    )


def build_design(L_q1_nH=L_Q1_NH, L_q2_nH=L_Q2_NH,
                 chip_size_x_mm=None, chip_size_y_mm=None,
                 chip_center_x_mm=None, chip_center_y_mm=None):
    """Build the two-transmon + interdigital capacitive coupler design.

    Q1 sits at ``pos_x = -D``, Q2 at ``pos_x = +D`` with ``D`` computed from the
    coupler geometry (:func:`_transmon_offset_um`) so each pad edge overlaps the
    coupler's lead stub by ``ABUT_OVERLAP_UM``. The coupler is rotated 90 deg so
    its comb-separation axis lies along global x: comb ``a`` abuts Q2 (+x), comb
    ``b`` abuts Q1 (-x).

    The chip is sized after the geometry is drawn, from the geometry's own
    bounding box * ``CHIP_MARGIN_FACTOR``, centred on that box. Any ``chip_*``
    kwarg overrides the corresponding auto-computed value.
    """
    design = designs.DesignPlanar({}, overwrite_enabled=True)

    design.variables['cpw_width'] = '20um'
    design.variables['cpw_gap'] = '12.25um'

    d_um = _transmon_offset_um()

    CircTransmon(design, Q1_NAME, options=_transmon_options(-d_um, L_q1_nH, CJ_Q1_FF))
    CircTransmon(design, Q2_NAME, options=_transmon_options(+d_um, L_q2_nH, CJ_Q2_FF))
    CapNInterdigitalRound(design, COUPLER_NAME, options=_coupler_options())

    # --- size the chip from the drawn geometry's own bounding box ---
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
    print(f"L_Q1_NH={L_Q1_NH}, L_Q2_NH={L_Q2_NH}, "
          f"transmon offset = +-{_transmon_offset_um():.1f}um")
    d = build_design()
    d.rebuild()
    print("Design built OK. Components:", list(d.components.keys()))
    print("Junction rows:", len(d.qgeometry.tables['junction']))
