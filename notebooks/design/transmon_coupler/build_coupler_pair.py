"""
Two CircTransmons joined by one CapNInterdigitalRound coupler, for the
coupler-geometry -> coupling-capacitance sweep.

vs build_chain.py this is a single coupler with the two pads pulled further
apart (a fixed NECK_UM of bare lead between each spine and its pad), so the
measured |C[Q1][Q2]| is dominated by the coupler itself and the direct
pad-to-pad term is small and (with the neck length held fixed) roughly
constant across the finger sweep.

The CapNInterdigitalRound API changed on 2026-09-01 (per-comb island_*_a/_b,
dog-leg lead_out_*/lead_across_*/end_*_turn, flare_overlap, corner_radius).
We use only the plain circular-island mode: island_radius + island_gap +
island_center_offset, no dog-leg. Pad centres come from coupler.island_centers().

Env / kwargs (all lengths um):
    FINGER_COUNT, FINGER_LENGTH_UM, FINGER_GAP_UM, FINGER_WIDTH_UM,
    FINGER_TIP_GAP_UM, SPINE_WIDTH_UM, NECK_UM, PAD_RADIUS_UM, PAD_GAP_UM,
    WITH_COUPLER (1/0 -- 0 draws only the two pads for a pad-pad reference)
"""
import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.qubits.circle_transmon import CircTransmon
from qiskit_metal.qlibrary.couplers.cap_n_interdigital_round import CapNInterdigitalRound

from build_design import _transmon_options, L_Q1_NH, CJ_Q1_FF, CHIP_SIZE_Z_UM

# ---- defaults = the 6-chain coupler ----
FINGER_COUNT       = 8
FINGER_LENGTH_UM   = 80.0
FINGER_GAP_UM      = 4.0
FINGER_WIDTH_UM    = 8.0
FINGER_TIP_GAP_UM  = 6.0
SPINE_WIDTH_UM     = 12.0
NECK_UM            = 60.0        # bare lead between spine outer edge and pad edge
PAD_RADIUS_UM      = 136.0
PAD_GAP_UM         = 20.0

CPL_FILLET_RADIUS_UM     = 6.0
CPL_CAP_GAP_GROUND_UM    = 3.0
CPL_ETCH_CORNER_RADIUS_UM = 16.0
CPL_LEAD_BLEND_RADIUS_UM = 40.0
CPL_ETCH_BLEND_RADIUS_UM = 10.0
CPL_FLARE_OVERLAP_UM     = 4.0
CPL_ROUND_RESOLUTION     = 32
Q1_NAME, Q2_NAME, COUPLER_NAME = "Q1", "Q2", "coupler"


def _island_center_offset(finger_length_um, finger_tip_gap_um, spine_width_um,
                          neck_um, pad_radius_um):
    """Pad-centre distance from the coupler centre such that a bare lead of
    length ``neck_um`` sits between the spine's outer edge and the pad edge.
    half_h = (fl + ftg)/2  (spine inner edge);  spine outer edge = half_h + sw.
    """
    half_h = 0.5 * (finger_length_um + finger_tip_gap_um)
    return half_h + spine_width_um + neck_um + pad_radius_um


def _coupler_options(pos_x_um, ico_um, finger_count, finger_length_um,
                     finger_gap_um, finger_width_um, finger_tip_gap_um,
                     spine_width_um, pad_radius_um, pad_gap_um):
    return Dict(
        pos_x=f"{pos_x_um}um", pos_y="0um", orientation="90", chip="main",
        finger_count=f"{finger_count}",
        finger_width=f"{finger_width_um}um",
        finger_length=f"{finger_length_um}um",
        finger_gap=f"{finger_gap_um}um",
        finger_tip_gap=f"{finger_tip_gap_um}um",
        spine_width=f"{spine_width_um}um",
        fillet_radius=f"{CPL_FILLET_RADIUS_UM}um",
        cap_gap_ground=f"{CPL_CAP_GAP_GROUND_UM}um",
        etch_corner_radius=f"{CPL_ETCH_CORNER_RADIUS_UM}um",
        lead_length="24um",
        lead_width=f"{spine_width_um}um",
        round_resolution=f"{CPL_ROUND_RESOLUTION}",
        island_radius=f"{pad_radius_um}um",
        island_gap=f"{pad_gap_um}um",
        island_center_offset=f"{ico_um}um",
        lead_flare_width=f"{3.0 * spine_width_um}um",
        lead_blend_radius=f"{CPL_LEAD_BLEND_RADIUS_UM}um",
        etch_blend_radius=f"{CPL_ETCH_BLEND_RADIUS_UM}um",
        flare_overlap=f"{CPL_FLARE_OVERLAP_UM}um",
    )


def build_coupler_pair(finger_count=FINGER_COUNT,
                       finger_length_um=FINGER_LENGTH_UM,
                       finger_gap_um=FINGER_GAP_UM,
                       finger_width_um=FINGER_WIDTH_UM,
                       finger_tip_gap_um=FINGER_TIP_GAP_UM,
                       spine_width_um=SPINE_WIDTH_UM,
                       neck_um=NECK_UM,
                       pad_radius_um=PAD_RADIUS_UM, pad_gap_um=PAD_GAP_UM,
                       with_coupler=True, chip_margin_um=700.0):
    design = designs.DesignPlanar({}, overwrite_enabled=True)
    design.variables["cpw_width"] = "20um"
    design.variables["cpw_gap"] = "12.25um"

    ico = _island_center_offset(finger_length_um, finger_tip_gap_um,
                                spine_width_um, neck_um, pad_radius_um)

    if with_coupler:
        cpl = CapNInterdigitalRound(
            design, COUPLER_NAME,
            options=_coupler_options(0.0, ico, finger_count, finger_length_um,
                                     finger_gap_um, finger_width_um,
                                     finger_tip_gap_um, spine_width_um,
                                     pad_radius_um, pad_gap_um))
        (xa, _), (xb, _) = cpl.island_centers()            # qiskit-metal p-units = mm
        xa, xb = xa * 1000.0, xb * 1000.0                  # -> um
    else:
        xa, xb = -ico, ico

    CircTransmon(design, Q1_NAME,
                 options=_transmon_options(xa, L_Q1_NH, CJ_Q1_FF))
    CircTransmon(design, Q2_NAME,
                 options=_transmon_options(xb, L_Q1_NH, CJ_Q1_FF))
    # override pad geometry (defaults in _transmon_options come from build_design)
    for q in (Q1_NAME, Q2_NAME):
        design.components[q].options.pad_radius = f"{pad_radius_um}um"
        design.components[q].options.pad_gap = f"{pad_gap_um}um"

    design.rebuild()

    # size the chip from the drawn bbox + a fixed margin
    allb = []
    for _, t in design.qgeometry.tables.items():
        if len(t):
            allb.append(t.total_bounds)
    allb = np.array(allb)
    minx, miny = allb[:, 0].min(), allb[:, 1].min()
    maxx, maxy = allb[:, 2].max(), allb[:, 3].max()
    m = chip_margin_um / 1000.0
    design.chips.main.size.size_x = f"{(maxx - minx) + 2 * m}mm"
    design.chips.main.size.size_y = f"{(maxy - miny) + 2 * m}mm"
    design.chips.main.size.size_z = f"{CHIP_SIZE_Z_UM}um"
    design.chips.main.size.center_x = f"{0.5 * (minx + maxx)}mm"
    design.chips.main.size.center_y = f"{0.5 * (miny + maxy)}mm"
    return design


if __name__ == "__main__":
    for wc in (True, False):
        d = build_coupler_pair(with_coupler=wc)
        d.rebuild()
        b = np.array([t.total_bounds for _, t in d.qgeometry.tables.items() if len(t)])
        print(f"with_coupler={wc}: components={list(d.components.keys())}  "
              f"chip={d.chips.main.size.size_x} x {d.chips.main.size.size_y}  "
              f"bbox_x=[{b[:,0].min():.3f},{b[:,2].max():.3f}]mm")
