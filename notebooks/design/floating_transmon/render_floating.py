"""Geometry check for the floating-transmon designs in build_design.py -- the
four layouts shown in the repo images. No meshing, no Palace."""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from render_geometry import render                       # noqa: E402
from build_design import (                               # noqa: E402
    single_transmon, coupled_pair, _design_bounds_mm,
    Q1_NAME, COUPLER_NAME,
)


def _islands(design, comp):
    t = design.qgeometry.tables["poly"]
    cid = design.components[comp].id
    L = unary_union(list(t[(t.component == cid) & (t.name == "island_L")]["geometry"]))
    R = unary_union(list(t[(t.component == cid) & (t.name == "island_R")]["geometry"]))
    return L, R


def _figure(design, fname, title, xlim=None, ylim=None):
    minx, miny, maxx, maxy = _design_bounds_mm(design)
    fig, ax = plt.subplots(figsize=(9, 8))
    render(design, ax)
    ax.set_xlim(*(xlim or (minx - 0.03, maxx + 0.03)))
    ax.set_ylim(*(ylim or (miny - 0.03, maxy + 0.03)))
    ax.set_title(title)
    fig.savefig(os.path.join(HERE, fname), dpi=170, bbox_inches="tight")
    plt.close(fig)


def _report(design, comp):
    L, R = _islands(design, comp)
    gap = L.distance(R) * 1e3
    print(f"  {comp}: islands valid={L.is_valid and R.is_valid}  disjoint={not L.intersects(R)}  "
          f"min L-R gap = {gap:.2f} um")


def main():
    d1 = single_transmon()
    _report(d1, "Q")
    _figure(d1, "floating_bare.png", "bare split + centre JJ")

    d2 = single_transmon(idc=True)
    _report(d2, "Q")
    _figure(d2, "floating_idc.png", "split + interdigitated shunt capacitor")

    d3 = single_transmon(idc=True, claw=True, charge=True)
    _report(d3, "Q")
    _figure(d3, "floating_full.png", "IDC + readout claw + charge line")

    dp = coupled_pair()
    _report(dp, Q1_NAME)
    t = dp.qgeometry.tables["poly"]
    sub = t["subtract"].astype(bool)
    metal = unary_union(list(t[~sub]["geometry"]))
    nets = [g for g in getattr(metal, "geoms", [metal]) if g.area * 1e6 > 300]
    worst = min((a.distance(b) for i, a in enumerate(nets) for b in nets[i + 1:]),
                default=float("nan"))
    print(f"  pair: {len(nets)} metal nets  min gap {worst * 1e3:.2f} um  "
          f"junctions {len(dp.qgeometry.tables['junction'])}")
    _figure(dp, "floating_pair.png", "two floating transmons + CapNInterdigitalRound coupler",
            xlim=(-0.55, 0.55), ylim=(-0.3, 0.3))

    D = dp.components[Q1_NAME].options["pos_x"]
    Dmm = -float(str(D).replace("um", "")) * 1e-3
    _figure(dp, "floating_junction.png", "Q1 -- coupler junction",
            xlim=(-Dmm - 0.02, -0.02), ylim=(-0.16, 0.16))

    print("saved floating_*.png")


if __name__ == "__main__":
    main()
