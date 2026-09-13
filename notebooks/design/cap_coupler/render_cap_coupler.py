"""Headless geometry check for the capacitive-coupler designs. Saves a PNG per
layout and prints net / ground / overlap counts. No meshing, no Palace."""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from shapely.ops import unary_union

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from render_geometry import render, sanity_checks            # noqa: E402
from build_design import coupled_pair, _design_bounds_mm, COUPLER_NAME  # noqa: E402


def _overlaps(design):
    t = design.qgeometry.tables["poly"]
    sub = t["subtract"].astype(bool)
    cid = design.components[COUPLER_NAME].id
    mc = unary_union(list(t[(~sub) & (t.component == cid)]["geometry"]))
    for nm in ("Q_gnd", "Q_float"):
        c = design.components[nm].id
        ov = unary_union(list(t[(~sub) & (t.component == c)]["geometry"])).intersection(mc).area * 1e6
        print(f"    coupler <-> {nm}: galvanic overlap {ov:.0f} um^2")


def _figure(design, fname, title):
    minx, miny, maxx, maxy = _design_bounds_mm(design)
    fig, ax = plt.subplots(figsize=(9, 10))
    render(design, ax)
    ax.set_xlim(minx - 0.03, maxx + 0.03)
    ax.set_ylim(miny - 0.03, maxy + 0.03)
    ax.set_title(title)
    fig.savefig(os.path.join(HERE, fname), dpi=160, bbox_inches="tight")
    plt.close(fig)


def main():
    for label, kw, fname in [
        ("direct", dict(shape="direct"), "cap_coupler_direct.png"),
        ("Z (dog-leg both)", dict(shape="Z"), "cap_coupler_Z.png"),
        ("U (dog-leg same side)", dict(shape="U"), "cap_coupler_U.png"),
        ("Z @ 30 deg", dict(shape="Z", orientation_deg=30), "cap_coupler_Z30.png"),
        ("Z, long legs", dict(shape="Z", lead_out_um=100, lead_across_um=80),
         "cap_coupler_long.png"),
    ]:
        d = coupled_pair(**kw)
        print(f"--- {label} ---")
        sanity_checks(d)
        _overlaps(d)
        _figure(d, fname, label)
    print("saved cap_coupler_*.png")


if __name__ == "__main__":
    main()
