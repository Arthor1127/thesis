"""A capacitive coupler between a grounded CircTransmon and a
CircTransmonFloating, with the optional dog-leg lead.

The coupler is ``CapNInterdigitalRound`` (env qlibrary/couplers) in
circular-island mode. By default each interdigital comb flares straight onto
its transmon pad (``shape='direct'``). With ``shape='Z'`` / ``'U'`` each comb
instead routes out on a rounded ``lead_out`` -> bend -> ``lead_across`` dog-leg
before flaring -- the fingers stay put and the two pads move off the direct
line between the qubits, so a second coupler (e.g. a ``SquidArrayCoupler``) can
share each island without the two interfering.

* ``shape='direct'`` -- comb flares straight onto the pad (the original coupler)
* ``shape='Z'``      -- combs dog-leg to opposite sides (pads diagonal)
* ``shape='U'``      -- combs dog-leg to the same side (pads on one side)

Layout is coupler-driven: the coupler is built first, then
``coupler.island_centers()`` says where the two transmons go.

* ``coupled_pair(...)``  -- grounded CircTransmon + CircTransmonFloating + coupler
* ``build_design()``     -- the default (Z) coupled pair
"""
import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.qubits.circle_transmon import CircTransmon
from qiskit_metal.qlibrary.qubits.circle_transmon_floating import CircTransmonFloating
from qiskit_metal.qlibrary.couplers.cap_n_interdigital_round import CapNInterdigitalRound
from qiskit_metal.qlibrary.couplers.squid_array_coupler import SquidArrayCoupler

# ---- transmon islands ----
GROUNDED_RADIUS_UM = 95.0
GROUNDED_GAP_UM = 10.0
FLOATING_RADIUS_UM = 110.0
FLOATING_GAP_UM = 12.0
FLOATING_ISLAND_GAP_UM = 6.0
GAP_CORNER_RADIUS_UM = 10.0

# ---- grounded-transmon junction ----
L_GND_NH = 11.0
CJ_GND_FF = 2.0

# ---- floating-transmon junction ----
L_FLOAT_NH = 12.0
CJ_FLOAT_FF = 2.0
FLOAT_JJ_WIDTH_UM = 6.0
FLOAT_JJ_SIM_GAP_UM = 3.0

# ---- CapNInterdigitalRound coupler ----
CPL_FINGER_COUNT = 8
CPL_FINGER_WIDTH_UM = 8.0
CPL_FINGER_LENGTH_UM = 70.0
CPL_FINGER_GAP_UM = 4.0
CPL_FINGER_TIP_GAP_UM = 6.0
CPL_SPINE_WIDTH_UM = 12.0
CPL_FILLET_RADIUS_UM = 4.0
CPL_CAP_GAP_GROUND_UM = 6.0
CPL_ETCH_CORNER_RADIUS_UM = 16.0
CPL_ETCH_BLEND_RADIUS_UM = 16.0
CPL_LEAD_LENGTH_UM = 20.0
CPL_LEAD_WIDTH_UM = 10.0
CPL_LEAD_FLARE_WIDTH_UM = 40.0
CPL_LEAD_BLEND_RADIUS_UM = 40.0
CPL_ROUND_RESOLUTION = 24
CPL_FLARE_OVERLAP_UM = 4.0
CPL_CORNER_RADIUS_UM = 18.0

# dog-leg leg lengths (shape != 'direct')
CPL_LEAD_OUT_UM = 60.0
CPL_LEAD_ACROSS_UM = 48.0
# direct-flare pad spacing (shape == 'direct'): half the pad-to-pad distance
ABUT_OVERLAP_UM = 4.0

# ---- names / chip ----
Q_GND_NAME, Q_FLOAT_NAME, COUPLER_NAME = 'Q_gnd', 'Q_float', 'cap_coupler'
CHIP_SIZE_Z_UM = 500
CHIP_MARGIN_FACTOR = 1.6
CPW_WIDTH = '20um'
CPW_GAP = '12um'

# per-shape (end_a_turn, end_b_turn); 'direct' -> no dog-leg
_SHAPE_TURNS = {'Z': (1, -1), 'U': (1, 1)}


def _direct_offset_um():
    """island_center_offset for shape='direct': the flare tip just reaches each
    pad rim (overlapping by ABUT_OVERLAP_UM)."""
    half_h = 0.5 * (CPL_FINGER_LENGTH_UM + CPL_FINGER_TIP_GAP_UM)
    return (GROUNDED_RADIUS_UM + half_h + CPL_SPINE_WIDTH_UM + CPL_LEAD_LENGTH_UM
            - ABUT_OVERLAP_UM)


def _design_bounds_mm(design):
    b = np.array([t.total_bounds for _, t in design.qgeometry.tables.items() if len(t)])
    return (b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max())


def _size_chip(design, margin=CHIP_MARGIN_FACTOR):
    minx, miny, maxx, maxy = _design_bounds_mm(design)
    cx, cy = 0.5 * (minx + maxx), 0.5 * (miny + maxy)
    design.chips.main.size.size_x = f'{(maxx - minx) * margin}mm'
    design.chips.main.size.size_y = f'{(maxy - miny) * margin}mm'
    design.chips.main.size.size_z = f'{CHIP_SIZE_Z_UM}um'
    design.chips.main.size.center_x = f'{cx}mm'
    design.chips.main.size.center_y = f'{cy}mm'


def coupler_options(orientation_deg=0.0, shape='Z', lead_out_um=CPL_LEAD_OUT_UM,
                    lead_across_um=CPL_LEAD_ACROSS_UM):
    o = Dict(
        pos_x='0um', pos_y='0um', orientation=f'{orientation_deg}', chip='main',
        finger_count=f'{CPL_FINGER_COUNT}',
        finger_width=f'{CPL_FINGER_WIDTH_UM}um',
        finger_length=f'{CPL_FINGER_LENGTH_UM}um',
        finger_gap=f'{CPL_FINGER_GAP_UM}um',
        finger_tip_gap=f'{CPL_FINGER_TIP_GAP_UM}um',
        spine_width=f'{CPL_SPINE_WIDTH_UM}um',
        fillet_radius=f'{CPL_FILLET_RADIUS_UM}um',
        cap_gap_ground=f'{CPL_CAP_GAP_GROUND_UM}um',
        etch_corner_radius=f'{CPL_ETCH_CORNER_RADIUS_UM}um',
        etch_blend_radius=f'{CPL_ETCH_BLEND_RADIUS_UM}um',
        lead_length=f'{CPL_LEAD_LENGTH_UM}um',
        lead_width=f'{CPL_LEAD_WIDTH_UM}um',
        lead_flare_width=f'{CPL_LEAD_FLARE_WIDTH_UM}um',
        lead_blend_radius=f'{CPL_LEAD_BLEND_RADIUS_UM}um',
        round_resolution=f'{CPL_ROUND_RESOLUTION}',
        island_radius=f'{GROUNDED_RADIUS_UM}um',      # shared fallback (direct mode)
        island_gap=f'{GROUNDED_GAP_UM}um',
        island_radius_a=f'{GROUNDED_RADIUS_UM}um',     # comb a -> grounded pad
        island_radius_b=f'{FLOATING_RADIUS_UM}um',     # comb b -> floating island
        island_gap_a=f'{GROUNDED_GAP_UM}um',
        island_gap_b=f'{FLOATING_GAP_UM}um',
        corner_radius=f'{CPL_CORNER_RADIUS_UM}um',
        flare_overlap=f'{CPL_FLARE_OVERLAP_UM}um',
    )
    if shape == 'direct':
        o['island_center_offset'] = f'{_direct_offset_um()}um'
    else:
        turn_a, turn_b = _SHAPE_TURNS[shape]
        o['island_center_offset'] = '0um'
        o['lead_out_a'] = f'{lead_out_um}um'
        o['lead_out_b'] = f'{lead_out_um}um'
        o['lead_across_a'] = f'{lead_across_um}um'
        o['lead_across_b'] = f'{lead_across_um}um'
        o['end_a_turn'] = f'{turn_a}'
        o['end_b_turn'] = f'{turn_b}'
    return o


def _new_design():
    d = designs.DesignPlanar({}, overwrite_enabled=True)
    d.variables['cpw_width'] = CPW_WIDTH
    d.variables['cpw_gap'] = CPW_GAP
    return d


def coupled_pair(shape='Z', orientation_deg=0.0, lead_out_um=CPL_LEAD_OUT_UM,
                 lead_across_um=CPL_LEAD_ACROSS_UM):
    """Grounded CircTransmon + CircTransmonFloating joined by the capacitive
    coupler. The transmons are placed where ``coupler.island_centers()`` says.

    The flare always approaches each pad along the component's local +y (comb a
    -> island A) / -y (comb b -> island B), so after the global rotation the
    contact angles are ``orientation + 90`` and ``orientation + 270``."""
    d = _new_design()
    CapNInterdigitalRound(d, COUPLER_NAME, options=coupler_options(
        orientation_deg, shape, lead_out_um, lead_across_um))
    cpl = d.components[COUPLER_NAME]

    (ax, ay), (bx, by) = cpl.island_centers()          # global, mm
    ang_a = orientation_deg + 90.0                      # flare hits pad A here
    ang_b = orientation_deg + 270.0                     # ... and pad B here

    CircTransmon(d, Q_GND_NAME, options=Dict(
        pos_x=f'{ax * 1e3}um', pos_y=f'{ay * 1e3}um', orientation='0', chip='main',
        cpw_width=CPW_WIDTH, cpw_gap=CPW_GAP,
        pad_radius=f'{GROUNDED_RADIUS_UM}um', pad_gap=f'{GROUNDED_GAP_UM}um',
        jj_options=Dict(jj_width='5um', jj_angle=f'{(ang_a + 180.0) % 360.0}',
                        L_j=f'{L_GND_NH}nH', C_j=f'{CJ_GND_FF}fF',
                        export_mask=False, jj_sim_gap='1um'),
    ))

    CircTransmonFloating(d, Q_FLOAT_NAME, options=Dict(
        pos_x=f'{bx * 1e3}um', pos_y=f'{by * 1e3}um', orientation=f'{ang_b}', chip='main',
        cpw_width=CPW_WIDTH, cpw_gap=CPW_GAP,
        pad_radius=f'{FLOATING_RADIUS_UM}um', pad_gap=f'{FLOATING_GAP_UM}um',
        gap_corner_radius=f'{GAP_CORNER_RADIUS_UM}um',
        island_gap=f'{FLOATING_ISLAND_GAP_UM}um',
        jj_options=Dict(jj_width=f'{FLOAT_JJ_WIDTH_UM}um', L_j=f'{L_FLOAT_NH}nH',
                        C_j=f'{CJ_FLOAT_FF}fF', export_mask=False,
                        jj_sim_gap=f'{FLOAT_JJ_SIM_GAP_UM}um'),
    ))

    d.rebuild()
    _size_chip(d)
    return d


def build_design():
    return coupled_pair(shape='Z')


# ==========================================================================
# Two FLOATING transmons coupled by BOTH couplers at once (parallel C + L)
# ==========================================================================
FLOAT_DIRECT_OVERLAP_UM = 4.0
SQ_N_SQUIDS = 5
SQ_L_J_NH = 10.0
SQ_CJ_FF = 1.0
SQ_JJ_WIDTH_UM = 2.0
SQ_LEAD_WIDTH_UM = 10.0
SQ_LEAD_OUT_UM = 55.0
SQ_CORNER_RADIUS_UM = 16.0
SQ_LEAD_BLEND_RADIUS_UM = 40.0
Q1_NAME, Q2_NAME = 'Q1', 'Q2'
CAP_NAME, SQUID_NAME = 'cap', 'squid'


def _float_direct_offset_um():
    half_h = 0.5 * (CPL_FINGER_LENGTH_UM + CPL_FINGER_TIP_GAP_UM)
    return (FLOATING_RADIUS_UM + half_h + CPL_SPINE_WIDTH_UM + CPL_LEAD_LENGTH_UM
            - FLOAT_DIRECT_OVERLAP_UM)


def _float_opts(pos_x_mm, pos_y_mm, orientation_deg):
    return Dict(
        pos_x=f'{pos_x_mm * 1e3}um', pos_y=f'{pos_y_mm * 1e3}um',
        orientation=f'{orientation_deg}', chip='main',
        cpw_width=CPW_WIDTH, cpw_gap=CPW_GAP,
        pad_radius=f'{FLOATING_RADIUS_UM}um', pad_gap=f'{FLOATING_GAP_UM}um',
        gap_corner_radius=f'{GAP_CORNER_RADIUS_UM}um',
        island_gap=f'{FLOATING_ISLAND_GAP_UM}um',
        jj_options=Dict(jj_width=f'{FLOAT_JJ_WIDTH_UM}um', L_j=f'{L_FLOAT_NH}nH',
                        C_j=f'{CJ_FLOAT_FF}fF', export_mask=False,
                        jj_sim_gap=f'{FLOAT_JJ_SIM_GAP_UM}um'),
    )


def _gnd_opts(pos_x_mm, pos_y_mm, jj_angle_deg=180.0, pad_radius_um=FLOATING_RADIUS_UM,
             pad_gap_um=FLOATING_GAP_UM):
    """Grounded CircTransmon: a single circular pad. ``jj_angle`` at 180 deg
    (west) keeps the JJ clear of the cap coupler (north / south) and the SQUID
    coupler (east)."""
    return Dict(
        pos_x=f'{pos_x_mm * 1e3}um', pos_y=f'{pos_y_mm * 1e3}um',
        orientation='0', chip='main', cpw_width=CPW_WIDTH, cpw_gap=CPW_GAP,
        pad_radius=f'{pad_radius_um}um', pad_gap=f'{pad_gap_um}um',
        jj_options=Dict(jj_width='5um', jj_angle=f'{jj_angle_deg}',
                        L_j=f'{L_GND_NH}nH', C_j=f'{CJ_GND_FF}fF',
                        export_mask=False, jj_sim_gap='1um'),
    )


def _squid_opts(pos_xy_mm, orientation_deg, lead_across_um, lead_out_um=SQ_LEAD_OUT_UM,
                n_squids=SQ_N_SQUIDS):
    return Dict(
        pos_x=f'{pos_xy_mm[0] * 1e3}um', pos_y=f'{pos_xy_mm[1] * 1e3}um',
        orientation=f'{orientation_deg}', chip='main',
        n_squids=f'{int(n_squids)}', L_j=f'{SQ_L_J_NH}nH', C_j=f'{SQ_CJ_FF}fF',
        jj_width=f'{SQ_JJ_WIDTH_UM}um', squid_gap='20um',
        plate1_width='4um', segment_a_length='13um', segment_a_width='4um',
        segment_b_length='13um', segment_b_width='4um', JJ_gap='3um',
        segment_c_width='4um', segment_d_length='6um', segment_d_width='9um',
        lead_width=f'{SQ_LEAD_WIDTH_UM}um',
        lead_out_a=f'{lead_out_um}um', lead_out_b=f'{lead_out_um}um',
        lead_across_a=f'{lead_across_um}um', lead_across_b=f'{lead_across_um}um',
        corner_radius=f'{SQ_CORNER_RADIUS_UM}um',
        end_a_turn='1', end_b_turn='1',                 # U: both islands one side
        island_radius_a=f'{FLOATING_RADIUS_UM}um', island_radius_b=f'{FLOATING_RADIUS_UM}um',
        island_gap_a=f'{FLOATING_GAP_UM}um', island_gap_b=f'{FLOATING_GAP_UM}um',
        flare_overlap='4um', lead_blend_radius=f'{SQ_LEAD_BLEND_RADIUS_UM}um',
        lead_gap='6um', ground_gap='4um', etch_blend_radius='14um',
        round_resolution=f'{CPL_ROUND_RESOLUTION}',
    )


def _solve_placement(a_loc, b_loc, target_a, target_b):
    """(pos_xy, orientation_deg) that maps a coupler's local island centres
    ``a_loc`` / ``b_loc`` onto ``target_a`` / ``target_b`` (a rigid rotate +
    translate). Also returns the spacing residual (mm) -- the fit is exact only
    if |b_loc - a_loc| == |target_b - target_a|."""
    a_loc, b_loc = np.asarray(a_loc), np.asarray(b_loc)
    ta, tb = np.asarray(target_a), np.asarray(target_b)
    dl, dt = b_loc - a_loc, tb - ta
    th = np.arctan2(dt[1], dt[0]) - np.arctan2(dl[1], dl[0])
    c, s = np.cos(th), np.sin(th)
    R = np.array([[c, -s], [s, c]])
    t = ta - R @ a_loc
    return (float(t[0]), float(t[1])), float(np.degrees(th)), float(abs(np.hypot(*dl) - np.hypot(*dt)))


def both_couplers(extra_gap_um=0.0, n_squids=SQ_N_SQUIDS, cap_lead_width_um=None,
                  grounded=False):
    """Two transmons coupled by BOTH couplers at once: a fixed interdigital
    capacitance and a flux-tunable SQUID-array inductance, in parallel. The cap
    coupler runs straight between the two transmons (vertically); the SQUID
    coupler dog-legs around the east side.

    ``grounded=False`` (default) -> two ``CircTransmonFloating``, each split
    turned +-45 deg so the cap coupler lands on ``island_L`` and the SQUID on
    ``island_R`` -- two independent coupling paths.
    ``grounded=True`` -> two ``CircTransmon``, a single pad each: both couplers
    merge onto the *same* pad from different sides (cap from N/S, SQUID from E),
    so C and L sit truly in parallel between the two pad nodes; the JJ is put at
    180 deg (west), clear of both.

    ``extra_gap_um`` pushes each transmon a further ``extra_gap_um`` from the
    centre; the cap coupler's lead stretches to follow (``cap_lead_width_um``
    widens it into a clean CPW-style trace for large gaps). Use it with a
    longer ``n_squids`` array so the SQUID chain still fits around the side.

    NOTE: the two couplers + two transmons enclose a loop, so the ground plane
    splits into an inner and an outer plate (``sanity_checks`` -> 2 pieces) --
    tie them together in the solver, exactly as in ``notebooks/design/
    floating_ring/`` section 5.
    """
    d = _new_design()

    # --- 1. cap coupler: vertical, direct flare, sets the transmon spacing ---
    co = coupler_options(orientation_deg=0.0, shape='direct')
    for k in ('island_radius', 'island_radius_a', 'island_radius_b'):
        co[k] = f'{FLOATING_RADIUS_UM}um'
    for k in ('island_gap', 'island_gap_a', 'island_gap_b'):
        co[k] = f'{FLOATING_GAP_UM}um'
    co['island_center_offset'] = f'{_float_direct_offset_um() + extra_gap_um}um'
    if cap_lead_width_um is not None:
        co['lead_width'] = f'{cap_lead_width_um}um'
        co['lead_flare_width'] = f'{max(cap_lead_width_um, CPL_LEAD_FLARE_WIDTH_UM)}um'
    CapNInterdigitalRound(d, CAP_NAME, options=co)
    (q1x, q1y), (q2x, q2y) = d.components[CAP_NAME].island_centers()
    spacing = np.hypot(q2x - q1x, q2y - q1y)

    # --- 2. the two transmons ---
    if grounded:
        CircTransmon(d, Q1_NAME, options=_gnd_opts(q1x, q1y))
        CircTransmon(d, Q2_NAME, options=_gnd_opts(q2x, q2y))
    else:
        CircTransmonFloating(d, Q1_NAME, options=_float_opts(q1x, q1y, -45.0))
        CircTransmonFloating(d, Q2_NAME, options=_float_opts(q2x, q2y, +45.0))

    # --- 3. SQUID coupler: U-shape, placement solved to hit the same two
    #        island centres; tune lead_across so its island spacing matches
    #        (clamped >= 0 -> if the chain alone is longer than the gap, widen
    #        the transmon spacing or shorten the array instead) ---
    la = 60.0
    for _ in range(8):
        probe = designs.DesignPlanar({}, overwrite_enabled=True)
        SquidArrayCoupler(probe, 'p', options=_squid_opts((0, 0), 0.0, la, n_squids=n_squids))
        a_loc, b_loc = probe.components['p'].island_centers()
        d_loc = np.hypot(b_loc[0] - a_loc[0], b_loc[1] - a_loc[1])
        if abs(d_loc - spacing) < 2e-3:
            break
        la += (spacing - d_loc) * 1e3 / 2.0            # each side takes half
        la = max(6.0, la)
    if abs(d_loc - spacing) > 5e-3:
        import warnings
        warnings.warn(f"SQUID island spacing {d_loc*1e3:.0f}um != transmon gap "
                      f"{spacing*1e3:.0f}um (n_squids={n_squids}); raise extra_gap_um.",
                      stacklevel=2)

    pos, ori, resid = _solve_placement(a_loc, b_loc, (q1x, q1y), (q2x, q2y))
    SquidArrayCoupler(d, SQUID_NAME, options=_squid_opts(pos, ori, la, n_squids=n_squids))

    d.rebuild()
    _size_chip(d)
    return d


if __name__ == "__main__":
    for shape in ('direct', 'Z', 'U'):
        d = coupled_pair(shape=shape)
        nj = len(d.qgeometry.tables['junction'])
        print(f"shape={shape:7s} components={list(d.components.keys())}  "
              f"junctions={nj}  chip={d.chips.main.size.size_x} x {d.chips.main.size.size_y}")
    d = both_couplers()
    print(f"both        components={list(d.components.keys())}  "
          f"junctions={len(d.qgeometry.tables['junction'])}  "
          f"chip={d.chips.main.size.size_x} x {d.chips.main.size.size_y}")
