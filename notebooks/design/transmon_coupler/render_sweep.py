"""Render every geometry in params_sweep.txt so the sweep can be eyeballed
before it runs. Full-chip panel per point + a zoom on the coupler.
Writes sweep_geometry.png and sweep_geometry_zoom.png."""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from build_coupler_pair import build_coupler_pair, _island_center_offset
import render_geometry as rg

try:
    from scipy.special import ellipk
    EPS0, EPS_R = 8.8541878128e-12, 11.45
    eps_eff = EPS0 * (1 + EPS_R) / 2

    def c_idc(N, L, w, s):
        eta = w / (w + s)
        k = np.sin(np.pi * eta / 2)
        return eps_eff * (N - 1) * (L * 1e-6) * ellipk(k ** 2) / ellipk(1 - k ** 2) * 1e15
except Exception:
    c_idc = lambda *a: float("nan")

pts = []
for ln in open(f"{HERE}/params_sweep.txt"):
    ln = ln.split()
    if len(ln) == 5:
        pts.append(tuple(float(x) for x in ln))

base = (8, 80, 4, 8, 6)


def label(p):
    fc, fl, fg, fw, ftg = p
    tags = []
    names = ["N", "L", "g", "w", "t"]
    for i, (v, b) in enumerate(zip(p, base)):
        tags.append(f"{names[i]}={v:g}" + ("*" if v != b else ""))
    return "  ".join(tags) + f"\n$C_{{IDC}}$={c_idc(fc, fl, fw, fg):.1f} fF"


ncol = 3
nrow = -(-len(pts) // ncol)

for zoom, fname in ((False, "sweep_geometry.png"), (True, "sweep_geometry_zoom.png")):
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 2.7 * nrow))
    for ax, p in zip(axes.flat, pts):
        fc, fl, fg, fw, ftg = p
        d = build_coupler_pair(finger_count=int(fc), finger_length_um=fl,
                               finger_gap_um=fg, finger_width_um=fw,
                               finger_tip_gap_um=ftg)
        d.rebuild()
        rg.render(d, ax)
        if zoom:
            ico = _island_center_offset(fl, ftg, 12.0, 60.0, 136.0) / 1000.0
            ax.set_xlim(-ico * 0.9, ico * 0.9)
            ax.set_ylim(-0.11, 0.11)
        else:
            ax.autoscale()
        ax.set_title(label(p), fontsize=8)
        ax.set_xlabel(""); ax.set_ylabel("")
        ax.tick_params(labelsize=6)
    for ax in axes.flat[len(pts):]:
        ax.axis("off")
    fig.suptitle("Coupler-geometry sweep — " + ("coupler zoom" if zoom else "full chip")
                 + f"  ({len(pts)} points; * = differs from baseline N8/L80/g4/w8/t6)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(f"{HERE}/{fname}", dpi=115)
    print("wrote", fname)
