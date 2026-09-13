"""
Render where AMR refined: load the per-cycle VOLUME meshes from a Palace AMR
run, slice at the metal plane, and plot local element size. Smaller elements
= where the error indicator sent refinement.

(The boundary paraview output keeps the seed surface mesh across all cycles,
so it does NOT show the refinement - the volume mesh is the only place it is
visible.)

Usage:  python render_amr_refinement.py <run_outputFiles_dir> <out_png>
"""
import os
import sys

import numpy as np
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pv.OFF_SCREEN = True

OF = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else "amr_refinement.png"
PV = os.path.join(OF, "paraview", "eigenmode")

# map paraview cycles to AMR stages via the .pvd timestep order
cycles = sorted(int(d[5:]) for d in os.listdir(PV) if d.startswith("Cycle"))
labels = {}
for k, c in enumerate(cycles):
    labels[c] = ("base mesh" if k == 0 else
                 f"after AMR iter {k}" if k < len(cycles) - 1 else "final")

# geometry window (from build_design bounds, a bit of margin)
XL, YL = (-1.45, -0.05), (-1.0, 0.20)

fig, axes = plt.subplots(1, len(cycles), figsize=(6.2 * len(cycles), 6.5))
if len(cycles) == 1:
    axes = [axes]
for ax, c in zip(axes, cycles):
    m = pv.read(os.path.join(PV, f"Cycle{c:06d}", "data.pvtu"))
    b = m.bounds
    zmid = 0.5 * (b[4] + b[5])
    sl = m.slice(normal="z", origin=(0, 0, zmid))
    sl = sl.compute_cell_sizes(length=False, area=True, volume=False)
    cc = sl.cell_centers()
    xy = cc.points
    esize_um = np.sqrt(np.abs(np.asarray(cc["Area"]))) * 1e6
    m_in = ((xy[:, 0] > XL[0] / 1e3) & (xy[:, 0] < XL[1] / 1e3) &
            (xy[:, 1] > YL[0] / 1e3) & (xy[:, 1] < YL[1] / 1e3))
    vmax = np.percentile(esize_um[m_in], 85) if m_in.any() else np.percentile(esize_um, 85)
    s = ax.scatter(xy[m_in, 0] * 1e3, xy[m_in, 1] * 1e3, c=esize_um[m_in],
                   s=3, cmap="turbo_r", vmin=0, vmax=vmax, linewidths=0)
    ax.set_title(f"{labels[c]}\n{m.n_cells:,} vol cells, slab {sl.n_cells:,} elems")
    ax.set_xlim(XL); ax.set_ylim(YL); ax.set_aspect("equal")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    plt.colorbar(s, ax=ax, label="local element size (um)", shrink=0.8)

plt.suptitle(f"AMR refinement map  ({OF})", fontsize=10)
plt.tight_layout()
plt.savefig(OUT, dpi=115)
print("wrote", OUT)
