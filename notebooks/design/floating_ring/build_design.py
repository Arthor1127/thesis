"""A ring of N floating (pocket) transmons, each capacitively coupled to its
two nearest neighbours by a ``CapNInterdigitalRound``.

The N transmons sit on the vertices of a regular N-gon. Every transmon is
turned tangent to the ring, so one half-island faces the next neighbour and the
other faces the previous one; a coupler bridges each edge, flaring tangentially
onto the facing arc of the two islands it joins. N = 6 is a hexagon.

Parametrised entry point::

    ring(n=6, side_um=None, pad_radius_um=120.0, idc=False)

* ``n``           -- number of transmons (>= 4; at n = 3 the bonds are 60 deg
                     off each island axis and the coupler flare would reach
                     over the split and short the two half-islands).
* ``side_um``     -- vertex-to-vertex spacing. ``None`` picks the value that
                     makes each coupler's tangent flare just overlap the two
                     rims (``ABUT_OVERLAP_UM``); pass a number to stretch or
                     squeeze the lattice (a warning fires if the flare then
                     no longer reaches).
* ``pad_radius_um`` -- half-island disc radius (transmon + coupler kept in sync).
* ``idc``         -- fill each transmon's split with the interdigitated shunt
                     capacitor.

The transmon / coupler shape constants are imported from the sibling
``floating_transmon`` design so the two folders never drift.
"""
import importlib.util
import os
import warnings

import numpy as np
from qiskit_metal import designs
from qiskit_metal.qlibrary.qubits.circle_transmon_floating import CircTransmonFloating
from qiskit_metal.qlibrary.couplers.cap_n_interdigital_round import CapNInterdigitalRound

# Reuse the sibling floating_transmon design (same file name, so load it by path
# under an alias rather than via ``import build_design``).
_FT = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "floating_transmon",
                                    "build_design.py"))
_spec = importlib.util.spec_from_file_location("_ft_build_design", _FT)
_ft = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ft)

floating_options = _ft.floating_options
coupler_options = _ft.coupler_options
_design_bounds_mm = _ft._design_bounds_mm
PAD_RADIUS_UM = _ft.PAD_RADIUS_UM
CHIP_SIZE_Z_UM = _ft.CHIP_SIZE_Z_UM
ABUT_OVERLAP_UM = _ft.ABUT_OVERLAP_UM
CPL_FINGER_LENGTH_UM = _ft.CPL_FINGER_LENGTH_UM
CPL_FINGER_TIP_GAP_UM = _ft.CPL_FINGER_TIP_GAP_UM
CPL_SPINE_WIDTH_UM = _ft.CPL_SPINE_WIDTH_UM
CPL_LEAD_LENGTH_UM = _ft.CPL_LEAD_LENGTH_UM

CPW_WIDTH = "20um"
CPW_GAP = "12um"
CHIP_MARGIN_FACTOR = 1.4
N_DEFAULT = 6


# --------------------------------------------------------------------------
def _coupler_reach_um():
    """Distance from the coupler centre to its flare tip -- everything beyond
    ``island_center_offset`` (half a spine plus the lead run plus half the
    finger channel)."""
    half_h = 0.5 * (CPL_FINGER_LENGTH_UM + CPL_FINGER_TIP_GAP_UM)
    return half_h + CPL_SPINE_WIDTH_UM + CPL_LEAD_LENGTH_UM


def matched_side_um(pad_radius_um=PAD_RADIUS_UM):
    """Vertex spacing at which each coupler's flare tip overlaps the two island
    rims by exactly ``ABUT_OVERLAP_UM``."""
    return 2.0 * (pad_radius_um + _coupler_reach_um() - ABUT_OVERLAP_UM)


def ring_vertices(n, side_um):
    """(x, y) in um of the N regular-polygon vertices, vertex 0 on +x."""
    r = side_um / (2.0 * np.sin(np.pi / n))
    return [(r * np.cos(2.0 * np.pi * k / n), r * np.sin(2.0 * np.pi * k / n))
            for k in range(n)], r


# --------------------------------------------------------------------------
def _new_design():
    d = designs.DesignPlanar({}, overwrite_enabled=True)
    d.variables["cpw_width"] = CPW_WIDTH
    d.variables["cpw_gap"] = CPW_GAP
    return d


def _size_chip(design, margin=CHIP_MARGIN_FACTOR):
    minx, miny, maxx, maxy = _design_bounds_mm(design)
    cx, cy = 0.5 * (minx + maxx), 0.5 * (miny + maxy)
    design.chips.main.size.size_x = f"{(maxx - minx) * margin}mm"
    design.chips.main.size.size_y = f"{(maxy - miny) * margin}mm"
    design.chips.main.size.size_z = f"{CHIP_SIZE_Z_UM}um"
    design.chips.main.size.center_x = f"{cx}mm"
    design.chips.main.size.center_y = f"{cy}mm"


def q_name(k):
    return f"Q{k}"


def coupler_name(k, n):
    return f"c{k}_{(k + 1) % n}"


# --------------------------------------------------------------------------
def ring(n=N_DEFAULT, side_um=None, pad_radius_um=PAD_RADIUS_UM, idc=False):
    """Build the N-transmon nearest-neighbour ring. See module docstring."""
    if n < 4:
        raise ValueError(
            "n must be >= 4: each bond is 180/n deg off the island axis, and at "
            "n = 3 that 60 deg puts the coupler flare over the split, shorting "
            "the two half-islands.")

    matched = matched_side_um(pad_radius_um)
    s = matched if side_um is None else float(side_um)
    if abs(s - matched) / matched > 0.05:
        warnings.warn(
            f"side_um={s:.1f} is {100 * (s - matched) / matched:+.0f}% off the "
            f"matched spacing {matched:.1f}um; the coupler flare may not meet "
            "the island rims (gap) or may bury too deep (short).", stacklevel=2)

    verts, r_circ = ring_vertices(n, s)
    off_deg = 180.0 / n                       # island axis <-> bond, per side
    d = _new_design()

    for k in range(n):
        psi = 360.0 * k / n
        o = floating_options(0.0, orientation=psi + 90.0, idc=idc)  # tangent
        o["pos_x"] = f"{verts[k][0]}um"
        o["pos_y"] = f"{verts[k][1]}um"
        if pad_radius_um != PAD_RADIUS_UM:
            o["pad_radius"] = f"{pad_radius_um}um"
        CircTransmonFloating(d, q_name(k), options=o)

    for k in range(n):
        a, b = verts[k], verts[(k + 1) % n]
        mid = (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]))
        beta = np.degrees(np.arctan2(b[1] - a[1], b[0] - a[0]))
        co = coupler_options(s / 2.0)          # island_center_offset = half a side
        co["pos_x"] = f"{mid[0]}um"
        co["pos_y"] = f"{mid[1]}um"
        co["orientation"] = f"{beta + 90.0}"
        if pad_radius_um != PAD_RADIUS_UM:
            co["island_radius"] = f"{pad_radius_um}um"
        CapNInterdigitalRound(d, coupler_name(k, n), options=co)

    d.rebuild()
    _size_chip(d)
    return d


def hexagon(idc=False):
    """The N = 6 ring (the headline layout)."""
    return ring(6, idc=idc)


def build_design():
    """Default design = the hexagon (matches the other design folders)."""
    return hexagon()


if __name__ == "__main__":
    for label, dd in [
        ("square (n=4)", ring(4)),
        ("hexagon (n=6)", ring(6)),
        ("octagon (n=8)", ring(8)),
        ("hexagon + shunt IDC", ring(6, idc=True)),
        ("hexagon, side +15%", ring(6, side_um=1.15 * matched_side_um())),
    ]:
        nq = sum(k.startswith("Q") for k in dd.components)
        nc = len(dd.components) - nq
        j = len(dd.qgeometry.tables["junction"])
        print(f"{label:24s} qubits={nq} couplers={nc} junctions={j}  "
              f"chip={dd.chips.main.size.size_x} x {dd.chips.main.size.size_y}")
