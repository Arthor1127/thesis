"""Floating (pocket) circular transmon designs.

A CircTransmonFloating is a circular disc split along its local y axis into a
LEFT and a RIGHT half-island, shunted at the centre by one Josephson junction.
Options here build the layouts shown in the repo images:

* single_transmon()                     -- bare split + centre JJ
* single_transmon(idc=True)             -- split filled with an interdigitated
                                           shunt capacitor (one comb per island)
* single_transmon(idc=..., claw=True,
                  charge=True)           -- + hanger-readout claw + XY drive line
* coupled_pair(idc=...)                  -- two floating transmons joined by a
                                           CapNInterdigitalRound coupler that
                                           flares tangentially onto one island

The chip is sized after the geometry is drawn, from its own bounding box.
"""
import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.qubits.circle_transmon_floating import CircTransmonFloating
from qiskit_metal.qlibrary.couplers.cap_n_interdigital_round import CapNInterdigitalRound

# ---- transmon disc ----
PAD_RADIUS_UM = 120.0
PAD_GAP_UM = 12.0                 # island rim -> ground pocket
GAP_CORNER_RADIUS_UM = 10.0       # fillet where the split edge meets the rim
ISLAND_GAP_PLAIN_UM = 6.0         # L/R straight-edge separation, bare split
ISLAND_GAP_IDC_UM = 30.0          # ... widened to fit the shunt-IDC fingers

# ---- junction (one JJ, always bridges the centre) ----
L_J_NH = 12.0
CJ_FF = 2.0
JJ_WIDTH_UM = 6.0
JJ_SIM_GAP_UM = 3.0

# ---- interdigitated shunt capacitor across the split ----
IDC_FINGER_COUNT = 0             # 0 = auto-fill the whole diameter
IDC_FINGER_WIDTH_UM = 6.0
IDC_FINGER_GAP_UM = 6.0
IDC_FINGER_TIP_GAP_UM = 6.0
IDC_FILLET_RADIUS_UM = 3.0
IDC_ROUND_RESOLUTION = 20

# ---- hanger-readout claw (capacitively coupled, on the RIGHT island) ----
CLAW_LOCATION_DEG = 0
CLAW_WIDTH_UM = 16.0
CLAW_GAP_UM = 8.0
CLAW_LENGTH_UM = 80.0
CLAW_GROUND_SPACING_UM = 6.0
CLAW_CONNECTOR_LENGTH_UM = 35.0

# ---- XY drive / charge line (on the LEFT island) ----
CHARGE_LOCATION_DEG = 180
CHARGE_PAD_RADIUS_UM = 14.0
CHARGE_PAD_GAP_UM = 8.0
CHARGE_GROUND_SPACING_UM = 6.0
CHARGE_LINE_LENGTH_UM = 130.0

# ---- CapNInterdigitalRound coupler (for coupled_pair) ----
CPL_FINGER_COUNT = 8
CPL_FINGER_WIDTH_UM = 8.0
CPL_FINGER_LENGTH_UM = 70.0
CPL_FINGER_GAP_UM = 4.0
CPL_FINGER_TIP_GAP_UM = 6.0
CPL_SPINE_WIDTH_UM = 12.0
CPL_FILLET_RADIUS_UM = 6.0
CPL_CAP_GAP_GROUND_UM = 6.0
CPL_ETCH_CORNER_RADIUS_UM = 16.0
CPL_ETCH_BLEND_RADIUS_UM = 18.0
CPL_LEAD_LENGTH_UM = 36.0
CPL_LEAD_WIDTH_UM = 12.0
CPL_LEAD_FLARE_WIDTH_UM = 44.0
CPL_LEAD_BLEND_RADIUS_UM = 40.0
CPL_ROUND_RESOLUTION = 24
ABUT_OVERLAP_UM = 4.0

# ---- names / chip ----
Q_NAME = 'Q'
Q1_NAME, Q2_NAME, COUPLER_NAME = 'Q1', 'Q2', 'coupler'
CHIP_SIZE_Z_UM = 500
CHIP_MARGIN_FACTOR = 2.5

CPW_WIDTH = '20um'
CPW_GAP = '12um'


def _transmon_offset_um():
    """|pos_x| of each transmon in coupled_pair so the coupler flare tip
    reaches each island rim (overlapping by ABUT_OVERLAP_UM)."""
    half_h = 0.5 * (CPL_FINGER_LENGTH_UM + CPL_FINGER_TIP_GAP_UM)
    return PAD_RADIUS_UM + (half_h + CPL_SPINE_WIDTH_UM + CPL_LEAD_LENGTH_UM) - ABUT_OVERLAP_UM


def _design_bounds_mm(design):
    """(minx, miny, maxx, maxy) in mm across every drawn qgeometry table."""
    b = np.array([t.total_bounds for _, t in design.qgeometry.tables.items() if len(t)])
    return (b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max())


def _size_chip(design, margin=CHIP_MARGIN_FACTOR):
    minx, miny, maxx, maxy = _design_bounds_mm(design)
    cx, cy = (minx + maxx) / 2, (miny + maxy) / 2
    design.chips.main.size.size_x = f'{(maxx - minx) * margin}mm'
    design.chips.main.size.size_y = f'{(maxy - miny) * margin}mm'
    design.chips.main.size.size_z = f'{CHIP_SIZE_Z_UM}um'
    design.chips.main.size.center_x = f'{cx}mm'
    design.chips.main.size.center_y = f'{cy}mm'


def floating_options(pos_x_um=0.0, orientation=0, idc=False, claw=False, charge=False):
    """Option Dict for one CircTransmonFloating."""
    o = Dict(
        pos_x=f'{pos_x_um}um', pos_y='0um', orientation=f'{orientation}', chip='main',
        cpw_width=CPW_WIDTH, cpw_gap=CPW_GAP,
        pad_radius=f'{PAD_RADIUS_UM}um', pad_gap=f'{PAD_GAP_UM}um',
        gap_corner_radius=f'{GAP_CORNER_RADIUS_UM}um',
        island_gap=f'{ISLAND_GAP_IDC_UM if idc else ISLAND_GAP_PLAIN_UM}um',
        jj_options=Dict(
            jj_width=f'{JJ_WIDTH_UM}um', L_j=f'{L_J_NH}nH', C_j=f'{CJ_FF}fF',
            export_mask=False, jj_sim_gap=f'{JJ_SIM_GAP_UM}um',
        ),
        shunt_idc=Dict(
            make_idc=idc,
            finger_count=f'{IDC_FINGER_COUNT}',
            finger_width=f'{IDC_FINGER_WIDTH_UM}um',
            finger_gap=f'{IDC_FINGER_GAP_UM}um',
            finger_tip_gap=f'{IDC_FINGER_TIP_GAP_UM}um',
            fillet_radius=f'{IDC_FILLET_RADIUS_UM}um',
            round_resolution=f'{IDC_ROUND_RESOLUTION}',
        ),
    )
    if claw:
        o['claw_a_options'] = Dict(
            make_claw=True, connector_location=f'{CLAW_LOCATION_DEG}',
            claw_width=f'{CLAW_WIDTH_UM}um', claw_gap=f'{CLAW_GAP_UM}um',
            claw_length=f'{CLAW_LENGTH_UM}um', ground_spacing=f'{CLAW_GROUND_SPACING_UM}um',
            connector_length=f'{CLAW_CONNECTOR_LENGTH_UM}um',
        )
    if charge:
        o['charge_line'] = Dict(
            make_charge_line=True, location=f'{CHARGE_LOCATION_DEG}',
            pad_radius=f'{CHARGE_PAD_RADIUS_UM}um', pad_gap=f'{CHARGE_PAD_GAP_UM}um',
            ground_spacing=f'{CHARGE_GROUND_SPACING_UM}um',
            line_length=f'{CHARGE_LINE_LENGTH_UM}um',
        )
    return o


def coupler_options(D_um):
    """Option Dict for the CapNInterdigitalRound coupler between the two
    transmons at +-D_um."""
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
        etch_blend_radius=f'{CPL_ETCH_BLEND_RADIUS_UM}um',
        lead_length=f'{CPL_LEAD_LENGTH_UM}um',
        lead_width=f'{CPL_LEAD_WIDTH_UM}um',
        lead_flare_width=f'{CPL_LEAD_FLARE_WIDTH_UM}um',
        lead_blend_radius=f'{CPL_LEAD_BLEND_RADIUS_UM}um',
        round_resolution=f'{CPL_ROUND_RESOLUTION}',
        island_radius=f'{PAD_RADIUS_UM}um',
        island_gap=f'{PAD_GAP_UM}um',
        island_center_offset=f'{D_um}um',
    )


def _new_design():
    d = designs.DesignPlanar({}, overwrite_enabled=True)
    d.variables['cpw_width'] = CPW_WIDTH
    d.variables['cpw_gap'] = CPW_GAP
    return d


def single_transmon(idc=False, claw=False, charge=False):
    """One CircTransmonFloating, centred at the origin."""
    d = _new_design()
    CircTransmonFloating(d, Q_NAME, options=floating_options(0.0, idc=idc, claw=claw, charge=charge))
    d.rebuild()
    _size_chip(d, margin=1.9 if (claw or charge) else 2.2)
    return d


def coupled_pair(idc=False, orientation_deg=0.0):
    """Two CircTransmonFloating at +-D joined by a CapNInterdigitalRound
    coupler flaring onto the facing island of each.

    ``orientation_deg`` rotates the whole 3-component assembly about the
    origin: the coupling axis then runs at that angle and each transmon is
    turned by the same amount, so its facing island still meets the coupler
    head-on. The layout is otherwise identical (same D, same coupler)."""
    d = _new_design()
    D = _transmon_offset_um()
    th = np.deg2rad(orientation_deg)
    c, s = np.cos(th), np.sin(th)

    for name, sign in ((Q1_NAME, -1.0), (Q2_NAME, +1.0)):
        o = floating_options(sign * D, orientation=orientation_deg, idc=idc)
        o['pos_x'] = f'{sign * D * c}um'
        o['pos_y'] = f'{sign * D * s}um'
        CircTransmonFloating(d, name, options=o)

    co = coupler_options(D)
    co['orientation'] = f'{90.0 + orientation_deg}'
    CapNInterdigitalRound(d, COUPLER_NAME, options=co)
    d.rebuild()
    _size_chip(d)
    return d


def build_design():
    """Default design = the coupled pair (matches the other design folders)."""
    return coupled_pair()


if __name__ == "__main__":
    for name, d in [("single", single_transmon()),
                    ("single+idc", single_transmon(idc=True)),
                    ("single+idc+claw+charge", single_transmon(idc=True, claw=True, charge=True)),
                    ("coupled_pair", coupled_pair()),
                    ("coupled_pair@35deg", coupled_pair(orientation_deg=35.0))]:
        j = len(d.qgeometry.tables['junction'])
        print(f"{name:26s} components={list(d.components.keys())}  junctions={j}  "
              f"chip={d.chips.main.size.size_x} x {d.chips.main.size.size_y}")
