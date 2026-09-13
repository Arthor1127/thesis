"""Geometry-only check: build the design and render it with matplotlib, NO
meshing and NO Palace. Produces a full-chip view plus a zoom on the coupler so
the finger rounding and the pad abutment can be reviewed cheaply."""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import build_design, _design_bounds_mm, COUPLER_NAME

HERE = os.path.dirname(os.path.abspath(__file__))


def _poly_patches(geom, **kw):
    """Yield matplotlib PathPatch(es) for a shapely (Multi)Polygon, holes included."""
    geoms = getattr(geom, "geoms", [geom])
    for g in geoms:
        if g.is_empty:
            continue
        verts, codes = [], []
        for ring in [g.exterior, *g.interiors]:
            xy = np.asarray(ring.coords)
            verts.extend(xy)
            codes.append(Path.MOVETO)
            codes.extend([Path.LINETO] * (len(xy) - 2))
            codes.append(Path.CLOSEPOLY)
        yield PathPatch(Path(verts, codes), **kw)


def render(design, ax):
    from shapely.ops import unary_union

    tables = design.qgeometry.tables
    poly = tables["poly"]
    sub = poly["subtract"].astype(bool)
    # dissolve all ground cutouts into one shape so overlapping subtract polygons
    # (e.g. the coupler moat re-drawing the transmon pockets) don't render as
    # spurious internal lines
    if sub.any():
        for p in _poly_patches(unary_union(list(poly[sub]["geometry"])),
                               fc="none", ec="crimson", hatch="///", lw=0.4, zorder=1):
            ax.add_patch(p)
    for _, row in poly[~sub].iterrows():
        for p in _poly_patches(row["geometry"],
                               fc="#3b7dd8", ec="#16375e", lw=0.6, zorder=2, alpha=0.9):
            ax.add_patch(p)
    if "junction" in tables and len(tables["junction"]):
        for _, row in tables["junction"].iterrows():
            xy = np.asarray(row["geometry"].coords)
            ax.plot(xy[:, 0], xy[:, 1], "-", color="limegreen", lw=2, zorder=3)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")


def sanity_checks(design):
    """Print net separation / ground connectivity so geometry bugs are loud."""
    from shapely.ops import unary_union
    from shapely.geometry import box as _box

    p = design.qgeometry.tables["poly"]
    sub = p["subtract"].astype(bool)

    metal = unary_union(list(p[~sub]["geometry"]))
    pieces = list(getattr(metal, "geoms", [metal]))
    big = sorted((g for g in pieces if g.area * 1e6 > 500), key=lambda g: -g.area)
    print(f"metal nets (>500 um2): {len(big)}")
    worst = min((a.distance(b) for i, a in enumerate(big) for b in big[i + 1:]),
               default=float("inf"))
    print(f"min gap between distinct metal nets: {worst * 1e3:.2f} um")

    sx = float(design.chips.main.size.size_x[:-2]); sy = float(design.chips.main.size.size_y[:-2])
    cx = float(design.chips.main.size.center_x[:-2]); cy = float(design.chips.main.size.center_y[:-2])
    chip = _box(cx - sx / 2, cy - sy / 2, cx + sx / 2, cy + sy / 2)
    ground = chip.difference(unary_union(list(p[sub]["geometry"])))
    gpieces = [g for g in getattr(ground, "geoms", [ground]) if g.area > 1e-7]
    print(f"ground pieces (>1e-7 mm2): {len(gpieces)}  [want 1 = fully connected]")


def main():
    design = build_design()
    design.rebuild()

    minx, miny, maxx, maxy = _design_bounds_mm(design)
    print(f"Design bounds (mm): X=[{minx:.4f},{maxx:.4f}]  Y=[{miny:.4f},{maxy:.4f}]")
    print(f"Chip (mm): {design.chips.main.size.size_x} x {design.chips.main.size.size_y}")
    print(f"Components: {list(design.components.keys())}")
    sanity_checks(design)

    fig, ax = plt.subplots(figsize=(11, 7))
    render(design, ax)
    pad = 0.05
    ax.set_xlim(minx - pad, maxx + pad)
    ax.set_ylim(miny - pad, maxy + pad)
    ax.set_title("transmon_coupler: full design")
    fig.savefig(os.path.join(HERE, "geometry_full.png"), dpi=160, bbox_inches="tight")

    # zoom on the coupler
    cg = design.components[COUPLER_NAME].qgeometry_bounds()  # (minx,miny,maxx,maxy) mm
    cx, cy = (cg[0] + cg[2]) / 2, (cg[1] + cg[3]) / 2
    half = max(cg[2] - cg[0], cg[3] - cg[1]) / 2 + 0.03
    fig, ax = plt.subplots(figsize=(9, 9))
    render(design, ax)
    ax.set_xlim(cx - half, cx + half)
    ax.set_ylim(cy - half, cy + half)
    ax.set_title("transmon_coupler: coupler zoom (rounded fingers + pad abutment)")
    fig.savefig(os.path.join(HERE, "geometry_coupler.png"), dpi=200, bbox_inches="tight")
    print("Saved geometry_full.png and geometry_coupler.png")


if __name__ == "__main__":
    main()
