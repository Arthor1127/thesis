"""
B_z(y) profile across the resonator hairpin, through a SQUID loop, for the
single-SQUID design and the array -- to see where the resonator flux is
screened.  Normalised to the resonator modal current so the two are comparable.

Usage: python compare_bz_profile.py <single_OF> <single_xc_um> <array_OF> <array_xc_um> [z_um]
"""
import sys
import numpy as np
import pandas as pd
import pyvista as pv
from scipy.special import ellipk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pv.OFF_SCREEN = True
MU0 = 4e-7 * np.pi
L0 = 1e-3
A_UM, B_UM = 10.0, 22.25
KK = float(ellipk((A_UM / B_UM) ** 2))
SQF = np.sqrt((((A_UM + B_UM) / 2) ** 2 - A_UM ** 2) * (B_UM ** 2 - ((A_UM + B_UM) / 2) ** 2)) * 1e-12
LEG = 50.0

sOF, sxc, aOF, axc = sys.argv[1], float(sys.argv[2]), sys.argv[3], float(sys.argv[4])
Z = float(sys.argv[5]) if len(sys.argv) > 5 else 2.5


def profile(OF, xc_um, label):
    vol = pv.read(f"{OF}/paraview/eigenmode/Cycle000001/data.pvtu")
    sl = vol.slice(normal="z", origin=(0, 0, Z * 1e-3))
    cc = sl.cell_centers().points
    slc = sl.point_data_to_cell_data()
    bz = np.asarray(slc["B_real"])[:, 2] + 1j * np.asarray(slc["B_imag"])[:, 2]
    # thin x-strip through the loop, bin in y
    m = np.abs(cc[:, 0] * 1e3 - xc_um) < 8.0
    y = cc[m, 1] * 1e3
    b = np.abs(bz[m])
    yb = np.arange(-80, 80, 1.5)
    prof = np.array([np.nanmean(b[(y >= yb[i]) & (y < yb[i + 1])]) if np.any((y >= yb[i]) & (y < yb[i + 1]))
                     else np.nan for i in range(len(yb) - 1)])
    yc = 0.5 * (yb[:-1] + yb[1:])
    # resonator current from the CPW gap field (both legs, both gaps)
    def bz_at(xx, yy, r=0.004):
        mm = (np.abs(cc[:, 0] - xx) < r) & (np.abs(cc[:, 1] - yy) < r)
        return complex(np.mean(bz[mm])) if mm.any() else np.nan
    vals = []
    for lg in (-LEG, LEG):
        for s in (-1, 1):
            v = bz_at(xc_um / 1e3, (lg + s * (A_UM + B_UM) / 2) / 1e3)
            if np.isfinite(v):
                vals.append(abs(v) * 4 * KK * SQF / (MU0 * B_UM * 1e-6))
    I_r = np.mean(vals)
    print(f"{label}: I_cpw(local) = {I_r:.3e} A   peak |Bz| in slot = {np.nanmax(prof[np.abs(yc)<20]):.3e} T")
    return yc, prof, I_r


ys, ps, Is = profile(sOF, sxc, "single SQUID")
ya, pa, Ia = profile(aOF, axc, "array loop")

fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
for a, norm, ttl in ((ax[0], 1.0, "raw |B_z| (mode-normalised, not comparable)"),
                     (ax[1], True, "|B_z| per unit resonator current  (T / A)")):
    sN = 1 / Is if norm is True else 1.0
    aN = 1 / Ia if norm is True else 1.0
    a.plot(ys, ps * sN, "-", color="#b9791f", lw=2, label="single SQUID (600 um loop)")
    a.plot(ya, pa * aN, "-", color="#1a887b", lw=2, label="array loop")
    a.axvspan(-8.75, 8.75, color="#b9791f", alpha=.10)
    a.axvspan(-12.5, 12.5, color="#1a887b", alpha=.10)
    for yy in (-50, 50):
        a.axvline(yy, color="#888", ls=":", lw=1)
    a.axvline(19.5, color="#1a887b", ls="--", lw=.8); a.axvline(-19.5, color="#1a887b", ls="--", lw=.8)
    a.set_xlabel("y across the hairpin (um)  -- resonator legs at +/- 50 (dotted)")
    a.set_ylabel("|B_z|" + (" / I_r  (T/A)" if norm is True else "  (a.u.)"))
    a.set_title(ttl, fontsize=10)
    a.set_xlim(-75, 75); a.grid(alpha=.3); a.legend(fontsize=8)
plt.suptitle(f"B_z across the hairpin through a SQUID loop   (z = {Z} um)\n"
             f"shaded = loop interior; array ground-pocket edge at +/-19.5 (dashed); single-SQUID pocket is +/-100 (no ground shown)")
plt.tight_layout()
plt.savefig("/tmp/bz_profile_compare.png", dpi=125)
print("saved /tmp/bz_profile_compare.png")

# integrated flux per unit x, per unit current  (the coupling density)
def flux_density(yc, prof, I_r, y_half):
    m = np.abs(yc) < y_half
    return np.trapz(prof[m], yc[m] * 1e-6) / I_r   # T*m / A  per unit x
print(f"\nflux per unit length per unit current (Wb/m / A):")
print(f"  single SQUID (|y|<8.75): {flux_density(ys, ps, Is, 8.75):.3e}")
print(f"  array loop   (|y|<12.5): {flux_density(ya, pa, Ia, 12.5):.3e}")
