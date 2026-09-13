"""Headless geometry check for the floating-transmon ring designs. Saves a PNG
per layout, prints net / ground / junction counts. No meshing, no Palace."""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from render_geometry import render, sanity_checks          # noqa: E402
from build_design import ring, hexagon, matched_side_um, _design_bounds_mm  # noqa: E402


def _figure(design, fname, title):
    minx, miny, maxx, maxy = _design_bounds_mm(design)
    fig, ax = plt.subplots(figsize=(9, 9))
    render(design, ax)
    ax.set_xlim(minx - 0.03, maxx + 0.03)
    ax.set_ylim(miny - 0.03, maxy + 0.03)
    ax.set_title(title)
    fig.savefig(os.path.join(HERE, fname), dpi=160, bbox_inches="tight")
    plt.close(fig)


def main():
    for label, design, fname in [
        ("hexagon (n=6)", hexagon(), "ring_hexagon.png"),
        ("square (n=4)", ring(4), "ring_n4.png"),
        ("octagon (n=8)", ring(8), "ring_n8.png"),
        ("hexagon + shunt IDC", ring(6, idc=True), "ring_hexagon_idc.png"),
        ("hexagon, side +20%", ring(6, side_um=1.2 * matched_side_um()),
         "ring_hexagon_wide.png"),
    ]:
        print(f"--- {label} ---")
        sanity_checks(design)
        _figure(design, fname, label)
    print("saved ring_*.png")


if __name__ == "__main__":
    main()
