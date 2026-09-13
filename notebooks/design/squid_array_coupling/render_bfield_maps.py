"""
Figures for the B-field mutual-inductance extraction:
 (a) |B| and B_z colormaps over one SQUID loop + the two resonator legs,
     at a few heights above the metal plane, with the flux-integration
     window drawn on;
 (b) where the field spikes (sharp PEC corners) and how little it changes
     the loop flux.

Usage: python render_bfield_maps.py <outputFiles_dir> [cycle] [loop] [outdir]
"""
import os
import sys
import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.patches import Rectangle

pv.OFF_SCREEN = True
ARRAY_X0_UM, PITCH_UM, LEG_Y_UM = -1290.0, 47.0, 50.0

OF = sys.argv[1]
CYC = int(sys.argv[2]) if len(sys.argv) > 2 else 1
LOOP = int(sys.argv[3]) if len(sys.argv) > 3 else 11
OUT = sys.argv[4] if len(sys.argv) > 4 else "/tmp"

f_r = pd.read_csv(f"{OF}/eig.csv").to_numpy()[CYC - 1, 1]
vol = pv.read(f"{OF}/paraview/eigenmode/Cycle{CYC:06d}/data.pvtu")
x0 = ARRAY_X0_UM + LOOP * PITCH_UM
xc = (x0 + 18) / 1e3                       # loop-interior centre, in mm
# flux window (um), same as extract_M_bfield.py
FW = dict(x=(x0 + 4, x0 + 37), y=(-12.5, 12.5))


def cut(z_um):
    sl = vol.slice(normal="z", origin=(0, 0, z_um * 1e-3))
    slc = sl.point_data_to_cell_data()
    cc = sl.cell_centers().points
    Bc = np.asarray(slc["B_real"]) + 1j * np.asarray(slc["B_imag"])
    win = (np.abs(cc[:, 0] - xc) < 0.06) & (np.abs(cc[:, 1]) < 0.075)
    return cc[win, 0] * 1e3, cc[win, 1] * 1e3, Bc[win]


# ---------- figure 1: |B| and B_z at z = 2.5 um ----------
x, y, B = cut(2.5)
Bmag = np.abs(B).sum(axis=1) if B.ndim > 1 else np.abs(B)
Bmag = np.linalg.norm(B, axis=1)
Bz = np.real(B[:, 2] * np.exp(-1j * np.angle(np.sum(B[:, 2]))))   # phase-aligned signed B_z

fig, ax = plt.subplots(1, 3, figsize=(20, 6.2))
sc = ax[0].scatter(x, y, c=Bmag, s=18, cmap="magma",
                   vmax=np.percentile(Bmag, 99.5), linewidths=0)
ax[0].set_title("|B|  (linear)")
plt.colorbar(sc, ax=ax[0], shrink=.82)
sc = ax[1].scatter(x, y, c=Bmag, s=18, cmap="magma",
                   norm=LogNorm(vmin=Bmag.max() / 300, vmax=Bmag.max()), linewidths=0)
ax[1].set_title("|B|  (log)")
plt.colorbar(sc, ax=ax[1], shrink=.82)
lim = np.percentile(np.abs(Bz), 98)
sc = ax[2].scatter(x, y, c=Bz, s=18, cmap="RdBu_r",
                   norm=TwoSlopeNorm(0, -lim, lim), linewidths=0)
ax[2].set_title("B_z  (signed, phase-aligned)")
plt.colorbar(sc, ax=ax[2], shrink=.82)
for a in ax:
    a.set_aspect("equal"); a.set_xlabel("x (um)"); a.set_ylabel("y (um)")
    a.axhline(LEG_Y_UM, color="c", lw=.7, ls=":"); a.axhline(-LEG_Y_UM, color="c", lw=.7, ls=":")
    a.add_patch(Rectangle((FW["x"][0], FW["y"][0]), FW["x"][1] - FW["x"][0],
                          FW["y"][1] - FW["y"][0], fill=False, ec="lime", lw=1.4, ls="--"))
plt.suptitle(f"resonator mode {f_r:.3f} GHz  -  SQUID loop {LOOP} + resonator legs (cyan dotted)  |  "
             f"z = 2.5 um above the metal plane  |  green dashed = flux-integration window")
plt.tight_layout()
p1 = f"{OUT}/bmap_loop{LOOP}_z2.5.png"
plt.savefig(p1, dpi=115); print("saved", p1)

# ---------- figure 2: height dependence + singularity ----------
fig, ax = plt.subplots(1, 4, figsize=(22, 5.6))
for a, z in zip(ax, [0.3, 1.0, 2.5, 6.0]):
    xx, yy, BB = cut(z)
    bm = np.linalg.norm(BB, axis=1)
    sc = a.scatter(xx, yy, c=bm, s=16, cmap="magma",
                   norm=LogNorm(vmin=bm.max() / 500, vmax=bm.max()), linewidths=0)
    a.set_title(f"z = {z} um   (|B| log,  max {bm.max():.1e})")
    a.set_aspect("equal"); a.set_xlabel("x (um)")
    a.add_patch(Rectangle((FW["x"][0], FW["y"][0]), FW["x"][1] - FW["x"][0],
                          FW["y"][1] - FW["y"][0], fill=False, ec="lime", lw=1.2, ls="--"))
    plt.colorbar(sc, ax=a, shrink=.8)
ax[0].set_ylabel("y (um)")
plt.suptitle(f"loop {LOOP}: |B| vs height above the plane -- the spikes sit on the metal edges / "
             f"junction-port corners (z<1um), OUTSIDE the green flux window")
plt.tight_layout()
p2 = f"{OUT}/bmap_loop{LOOP}_height.png"
plt.savefig(p2, dpi=115); print("saved", p2)

# ---------- quantify: flux with vs without the brightest cells, per height ----------
print(f"\n{'z(um)':>6} {'Phi (full)':>13} {'Phi (drop top 1%)':>18} {'max/median in window':>22}")
for z in [1.0, 2.0, 2.5, 3.0, 4.0, 6.0]:
    sl = vol.slice(normal="z", origin=(0, 0, z * 1e-3)).compute_cell_sizes(area=True)
    slc = sl.point_data_to_cell_data()
    cc = sl.cell_centers().points
    bz = np.asarray(slc["B_real"])[:, 2] + 1j * np.asarray(slc["B_imag"])[:, 2]
    A = np.asarray(sl["Area"])
    m = ((cc[:, 0] * 1e3 > FW["x"][0]) & (cc[:, 0] * 1e3 < FW["x"][1])
         & (np.abs(cc[:, 1] * 1e3) < 12.5))
    full = np.real(np.sum(bz[m] * A[m]))
    thr = np.percentile(np.abs(bz[m]), 99)
    kp = m & (np.abs(bz) <= thr)
    clip = np.real(np.sum(bz[kp] * A[kp]))
    print(f"{z:>6.1f} {full:>13.3e} {clip:>18.3e} {np.abs(bz[m]).max()/np.median(np.abs(bz[m])):>22.1f}")
