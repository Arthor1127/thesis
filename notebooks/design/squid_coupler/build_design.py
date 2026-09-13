"""A tunable SQUID-array coupler between a grounded CircTransmon and a
CircTransmonFloating.

The coupler (``SquidArrayCoupler``, in the env's qlibrary/couplers) is a
straight series chain of ``n_squids`` SQUID loops -- a flux-tunable
superinductor -- with a rounded dog-leg lead on each end. Each lead leaves its
transmon island tangentially (the same sweeping flare as the capacitive
``CapNInterdigitalRound``), runs *outward* from the island (``lead_out``), turns
through a rounded bend, runs back in (``lead_across``), and meets the array.
Extending outward first leaves room on each island for a second coupler.

``shape='Z'`` puts the two transmons on opposite sides of the array (a dog-leg
each way); ``shape='U'`` puts them on the same side (both leads dip to the
array and back).

Layout is coupler-driven: the coupler is built first, then
``coupler.island_centers()`` says where the two transmons must sit for the
flares to land tangent on their rims.

* ``coupled_pair(...)``  -- grounded CircTransmon  +  CircTransmonFloating  +  coupler
* ``build_design()``     -- the default coupled pair
"""
import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.qubits.circle_transmon import CircTransmon
from qiskit_metal.qlibrary.qubits.circle_transmon_floating import CircTransmonFloating
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
GND_JJ_ANGLE_DEG = 180          # points away from the coupler (coupler hits north)

# ---- floating-transmon junction ----
L_FLOAT_NH = 12.0
CJ_FLOAT_FF = 2.0
FLOAT_JJ_WIDTH_UM = 6.0
FLOAT_JJ_SIM_GAP_UM = 3.0

# ---- SQUID-array coupler ----
CPL_N_SQUIDS = 6
CPL_L_J_NH = 10.0               # per junction (2 * n_squids junctions total)
CPL_CJ_FF = 1.0
CPL_JJ_WIDTH_UM = 2.0
CPL_SQUID_GAP_UM = 22.0
CPL_ARM_WIDTH_UM = 4.0
CPL_ARM_SEMI_LEN_UM = 14.0
CPL_JJ_GAP_UM = 3.0
CPL_SEG_D_LEN_UM = 6.0
CPL_SEG_D_WID_UM = 9.0
CPL_LEAD_WIDTH_UM = 10.0
CPL_LEAD_OUT_UM = 70.0
CPL_LEAD_ACROSS_UM = 45.0
CPL_CORNER_RADIUS_UM = 18.0
CPL_FLARE_OVERLAP_UM = 4.0
CPL_LEAD_BLEND_RADIUS_UM = 44.0
CPL_LEAD_GAP_UM = 6.0
CPL_GROUND_GAP_UM = 4.0
CPL_ETCH_BLEND_RADIUS_UM = 14.0
CPL_ROUND_RESOLUTION = 24

# ---- names / chip ----
Q_GND_NAME, Q_FLOAT_NAME, COUPLER_NAME = 'Q_gnd', 'Q_float', 'squid_coupler'
CHIP_SIZE_Z_UM = 500
CHIP_MARGIN_FACTOR = 1.6
CPW_WIDTH = '20um'
CPW_GAP = '12um'

# per-end bend direction for each shape: (end_a_turn, end_b_turn)
_SHAPE_TURNS = {'Z': (-1, +1), 'U': (+1, +1), 'U_down': (-1, -1)}


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


def coupler_options(orientation_deg=0.0, shape='Z', n_squids=CPL_N_SQUIDS,
                    lead_out_um=CPL_LEAD_OUT_UM, lead_across_um=CPL_LEAD_ACROSS_UM):
    turn_a, turn_b = _SHAPE_TURNS[shape]
    return Dict(
        pos_x='0um', pos_y='0um', orientation=f'{orientation_deg}', chip='main',
        n_squids=f'{int(n_squids)}',
        L_j=f'{CPL_L_J_NH}nH', C_j=f'{CPL_CJ_FF}fF', jj_width=f'{CPL_JJ_WIDTH_UM}um',
        squid_gap=f'{CPL_SQUID_GAP_UM}um',
        plate1_width=f'{CPL_ARM_WIDTH_UM}um',
        segment_a_length=f'{CPL_ARM_SEMI_LEN_UM}um', segment_a_width=f'{CPL_ARM_WIDTH_UM}um',
        segment_b_length=f'{CPL_ARM_SEMI_LEN_UM}um', segment_b_width=f'{CPL_ARM_WIDTH_UM}um',
        JJ_gap=f'{CPL_JJ_GAP_UM}um', segment_c_width=f'{CPL_ARM_WIDTH_UM}um',
        segment_d_length=f'{CPL_SEG_D_LEN_UM}um', segment_d_width=f'{CPL_SEG_D_WID_UM}um',
        lead_width=f'{CPL_LEAD_WIDTH_UM}um',
        lead_out_a=f'{lead_out_um}um', lead_out_b=f'{lead_out_um}um',
        lead_across_a=f'{lead_across_um}um', lead_across_b=f'{lead_across_um}um',
        corner_radius=f'{CPL_CORNER_RADIUS_UM}um',
        end_a_turn=f'{turn_a}', end_b_turn=f'{turn_b}',
        island_radius_a=f'{GROUNDED_RADIUS_UM}um', island_gap_a=f'{GROUNDED_GAP_UM}um',
        island_radius_b=f'{FLOATING_RADIUS_UM}um', island_gap_b=f'{FLOATING_GAP_UM}um',
        flare_overlap=f'{CPL_FLARE_OVERLAP_UM}um',
        lead_blend_radius=f'{CPL_LEAD_BLEND_RADIUS_UM}um',
        lead_gap=f'{CPL_LEAD_GAP_UM}um', ground_gap=f'{CPL_GROUND_GAP_UM}um',
        etch_blend_radius=f'{CPL_ETCH_BLEND_RADIUS_UM}um',
        round_resolution=f'{CPL_ROUND_RESOLUTION}',
    )


def _contact_angle_deg(island_xy, bend_xy):
    """Angle (deg) from the island centre toward the point the lead comes from
    (the bend) -- i.e. where the flare lands on the rim."""
    return np.degrees(np.arctan2(bend_xy[1] - island_xy[1], bend_xy[0] - island_xy[0]))


def _new_design():
    d = designs.DesignPlanar({}, overwrite_enabled=True)
    d.variables['cpw_width'] = CPW_WIDTH
    d.variables['cpw_gap'] = CPW_GAP
    return d


def coupled_pair(shape='Z', orientation_deg=0.0, n_squids=CPL_N_SQUIDS,
                 lead_out_um=CPL_LEAD_OUT_UM, lead_across_um=CPL_LEAD_ACROSS_UM):
    """Grounded CircTransmon + CircTransmonFloating joined by the SQUID-array
    coupler. The transmons are placed where ``coupler.island_centers()`` says."""
    d = _new_design()
    cpl = SquidArrayCoupler(d, COUPLER_NAME, options=coupler_options(
        orientation_deg, shape, n_squids, lead_out_um, lead_across_um))

    (ax, ay), (bx, by) = cpl.island_centers()          # global, mm
    # each end's bend point, global mm, to get the flare contact angle
    p = cpl.p
    th = np.deg2rad(float(p.orientation))
    c, s = np.cos(th), np.sin(th)
    _, _, _, _, xl, xr, ch = cpl._build_chain()
    turn_a, turn_b = _SHAPE_TURNS[shape]
    xbend_a = xl + (-1) * p.lead_across_a
    xbend_b = xr + (+1) * p.lead_across_b
    bend_a = (c * xbend_a - s * 0.0 + float(p.pos_x), s * xbend_a + c * 0.0 + float(p.pos_y))
    bend_b = (c * xbend_b - s * 0.0 + float(p.pos_x), s * xbend_b + c * 0.0 + float(p.pos_y))

    ang_a = _contact_angle_deg((ax, ay), bend_a)       # coupler hits grounded pad here
    ang_b = _contact_angle_deg((bx, by), bend_b)       # ... and the floating island here

    # grounded transmon: circular pad, JJ pointed away from the contact
    CircTransmon(d, Q_GND_NAME, options=Dict(
        pos_x=f'{ax * 1e3}um', pos_y=f'{ay * 1e3}um', orientation='0', chip='main',
        cpw_width=CPW_WIDTH, cpw_gap=CPW_GAP,
        pad_radius=f'{GROUNDED_RADIUS_UM}um', pad_gap=f'{GROUNDED_GAP_UM}um',
        jj_options=Dict(jj_width='5um', jj_angle=f'{(ang_a + 180.0) % 360.0}',
                        L_j=f'{L_GND_NH}nH', C_j=f'{CJ_GND_FF}fF',
                        export_mask=False, jj_sim_gap='1um'),
    ))

    # floating transmon: one half-island's arc must face the contact -> set its
    # orientation so island_R (local +x, faces `orientation`) points at ang_b
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
    return coupled_pair()


if __name__ == "__main__":
    for shape in ('Z', 'U'):
        d = coupled_pair(shape=shape)
        nj = len(d.qgeometry.tables['junction'])
        print(f"shape={shape:3s} components={list(d.components.keys())}  "
              f"junction rows={nj}  chip={d.chips.main.size.size_x} x {d.chips.main.size.size_y}")
