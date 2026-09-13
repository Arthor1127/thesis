"""
Linear chain of N grounded circular transmons (CircTransmon), each adjacent
pair joined by one fixed interdigital capacitive coupler (CapNInterdigitalRound)
built exactly as in build_design.py: each comb is a galvanic extension of one
transmon island, the two combs of a coupler never touch.

Geometry is the 2-transmon build_design.py repeated on a uniform pitch:

    pitch  P = 2 * D,   D = _transmon_offset_um()      (build_design.py)
    transmon k   at  x = (k - (N-1)/2) * P            k = 0 .. N-1
    coupler   j   at  x = midpoint(x_j, x_{j+1})       j = 0 .. N-2

Every coupler sees an identical local geometry (two islands at +-D, same
options), so _coupler_options() is reused verbatim apart from pos_x.

Chip is sized from the drawn bounding box after everything is placed.
"""
from qiskit_metal import designs
from qiskit_metal.qlibrary.qubits.circle_transmon import CircTransmon
from qiskit_metal.qlibrary.couplers.cap_n_interdigital_round import CapNInterdigitalRound

from build_design import (
    L_Q1_NH, CJ_Q1_FF,
    CHIP_SIZE_Z_UM, CHIP_MARGIN_FACTOR,
    _design_bounds_mm, _transmon_offset_um,
    _transmon_options, _coupler_options,
)

N_TRANSMONS = 6


def _coupler_options_at(pos_x_um):
    """build_design._coupler_options() with pos_x overridden. The coupler's
    island_center_offset is already _transmon_offset_um() = the chain half-pitch,
    so the same options describe every coupler in the chain."""
    opts = _coupler_options()
    opts["pos_x"] = f"{pos_x_um}um"
    return opts


def transmon_name(k):
    return f"Q{k + 1}"


def coupler_name(j):
    return f"coupler_{j + 1}"


def build_chain(n_transmons=N_TRANSMONS, L_nH=L_Q1_NH,
                chip_size_x_mm=None, chip_size_y_mm=None,
                chip_center_x_mm=None, chip_center_y_mm=None):
    """Build the N-transmon capacitive-coupler chain.

    Returns the qiskit-metal DesignPlanar. Components are named Q1..QN
    (west -> east) and coupler_1..coupler_{N-1} (coupler_j sits between Qj
    and Q{j+1}).
    """
    design = designs.DesignPlanar({}, overwrite_enabled=True)
    design.variables["cpw_width"] = "20um"
    design.variables["cpw_gap"] = "12.25um"

    D = _transmon_offset_um()          # half-pitch (um)
    P = 2.0 * D                        # transmon pitch (um)
    x0 = -(n_transmons - 1) / 2.0 * P  # x of Q1, so the chain is centred on 0

    x_transmon = [x0 + k * P for k in range(n_transmons)]
    for k, xk in enumerate(x_transmon):
        CircTransmon(design, transmon_name(k),
                     options=_transmon_options(xk, L_nH, CJ_Q1_FF))

    for j in range(n_transmons - 1):
        x_mid = 0.5 * (x_transmon[j] + x_transmon[j + 1])
        CapNInterdigitalRound(design, coupler_name(j),
                              options=_coupler_options_at(x_mid))

    # --- size the chip from the drawn geometry's own bounding box ---
    minx, miny, maxx, maxy = _design_bounds_mm(design)
    bcx, bcy = (minx + maxx) / 2, (miny + maxy) / 2
    bsx, bsy = (maxx - minx), (maxy - miny)

    sx = bsx * CHIP_MARGIN_FACTOR if chip_size_x_mm is None else chip_size_x_mm
    sy = bsy * CHIP_MARGIN_FACTOR if chip_size_y_mm is None else chip_size_y_mm
    # the chain is long and thin: don't let the y-margin collapse to a sliver
    sy = max(sy, bsy + 1.5)
    cx = bcx if chip_center_x_mm is None else chip_center_x_mm
    cy = bcy if chip_center_y_mm is None else chip_center_y_mm

    design.chips.main.size.size_x = f"{sx}mm"
    design.chips.main.size.size_y = f"{sy}mm"
    design.chips.main.size.size_z = f"{CHIP_SIZE_Z_UM}um"
    design.chips.main.size.center_x = f"{cx}mm"
    design.chips.main.size.center_y = f"{cy}mm"
    return design


if __name__ == "__main__":
    d = build_chain()
    d.rebuild()
    D = _transmon_offset_um()
    print(f"N = {N_TRANSMONS} transmons, pitch = {2*D:.1f} um")
    print("components:", list(d.components.keys()))
    print("junction rows:", len(d.qgeometry.tables["junction"]))
    b = _design_bounds_mm(d)
    print(f"design bounds (mm): X=[{b[0]:.3f}, {b[2]:.3f}]  Y=[{b[1]:.3f}, {b[3]:.3f}]")
    print(f"chip (mm): {d.chips.main.size.size_x} x {d.chips.main.size.size_y}")
