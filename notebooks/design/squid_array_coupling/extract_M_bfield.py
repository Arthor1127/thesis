"""
Direct B-field extraction of the per-loop mutual inductance M_k and coupling
c_k, bypassing Palace's port-EPR (which floors out for the interior loops of
this weak coupling).

Method (user, 2026-08-29):

  1. Flux through loop k:  Phi_k = integral of B_z over the loop interior, on
     a horizontal cut a few um ABOVE the metal plane (on z=0 itself B.n = 0
     by the PEC BC, so the flux must be read just above it).

  2. Resonator modal current, from its OWN CPW gap field using the exact
     zero-thickness-PEC CPW result (Lahl & Woerdenweber, APL 81 505 (2002);
     Bothner et al. arXiv:1202.5963), lambda_eff -> 0:

        center |x|<a, gaps a<|x|<b, k = a/b, K = complete elliptic integral 1st kind
        B_z,gap(x) = mu0 I b / ( 4 K(k) sqrt( (x^2 - a^2)(b^2 - x^2) ) ),   a<|x|<b

     Invert at a sampled gap point x_g to get the local current
        I(x_along) = B_z,gap(x_g) * 4 K(k) sqrt((x_g^2-a^2)(b^2-x_g^2)) / (mu0 b)
     (cross-checked against the crude twin-lead value.)

  3. M_k = Phi_k / I_r , with I_r the modal current amplitude (peak of the
     standing wave). c_k = pi M_k I_r^zpf / Phi0, I_r^zpf = sqrt(hbar w_r/2 L_r),
     L_r = 2 U_m_tot / I_r^2.  Palace's U_m density is (|Br|^2+|Bi|^2)/(2 mu0)
     (verified: [(|Br|^2+|Bi|^2)/(4 mu0)] / U_m = 0.5 exactly), so the volume
     integral U_m_tot == domain-E.csv E_mag is the PEAK inductive energy
     (1/2) L_r I_r^2, not the time average.  I_r cancels:

        c_k = pi * Phi_k * sqrt( hbar w_r / (4 U_m_tot) ) / Phi0

     so the coupling needs only Phi_k and U_m_tot (volume integral of the U_m
     density). The CPW formula is used only to report M_k in henries and to
     sanity-check.  Both Phi_k and U_m_tot carry the arbitrary eigenmode
     normalisation, which cancels in the ratio -> c_k, M_k, L_r are all
     independent of the absolute field scale (verified over 9 orders).

UNITS (important -- fixed 2026-08-29 after a factor-1e6/1e9 bug):
  the Palace mesh has L0 = 1e-3, i.e. mesh coordinates are in MILLIMETRES.
  pyvista compute_cell_sizes therefore returns AREA in mm^2 and VOLUME in
  mm^3.  The exported B (Tesla) and U_m (J/m^3) fields are in SI.  So:
      area_SI  [m^2] = area_mesh [mm^2] * (L0)^2 = area_mesh * 1e-6
      vol_SI   [m^3] = vol_mesh  [mm^3] * (L0)^3 = vol_mesh  * 1e-9
  Cross-check: sum(U_m * dV_SI) must equal E_mag in outputFiles/domain-E.csv.

Usage:  python extract_M_bfield.py <outputFiles_dir> [cycle] [z_offset_um]
Writes: <outputFiles_dir>/../M_bfield_cyc{C}_z{Z}.npz   (raw arrays)
        /tmp/M_bfield.png                                (bar charts)
"""
import os
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
HBAR = 1.054571817e-34
PHI0 = 2.067833848e-15

L0 = 1.0e-3                       # Palace Model.L0: 1 mesh unit = 1 mm
AREA_TO_SI = L0 ** 2              # mm^2 -> m^2   (1e-6)
VOL_TO_SI = L0 ** 3               # mm^3 -> m^3   (1e-9)

# geometry (um)
ARRAY_X0_UM = -1290.0
PITCH_UM = 47.0
CPW_W_UM = 20.0
CPW_GAP_UM = 12.25
HAIRPIN_LEG_Y_UM = 50.0          # each leg centre is +/- this in y
A_UM = CPW_W_UM / 2              # 10
B_UM = A_UM + CPW_GAP_UM        # 22.25
KMOD = A_UM / B_UM
KK = float(ellipk(KMOD ** 2))    # scipy takes m = k^2

OF = sys.argv[1]
CYC = int(sys.argv[2]) if len(sys.argv) > 2 else 1
Z_OFF_UM = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0

f_r = pd.read_csv(f"{OF}/eig.csv").to_numpy()[CYC - 1, 1] * 1e9
w_r = 2 * np.pi * f_r
P = pd.read_csv(f"{OF}/port-EPR.csv").to_numpy()[:, 1:]
n_loops = P.shape[1] // 2
ploop_csv = P[CYC - 1][0::2] + P[CYC - 1][1::2]
print(f"cycle {CYC}: f_r = {f_r/1e9:.4f} GHz, {n_loops} loops, K(k={KMOD:.3f}) = {KK:.4f}")

vol = pv.read(f"{OF}/paraview/eigenmode/Cycle{CYC:06d}/data.pvtu").compute_cell_sizes(volume=True)
vc = vol.point_data_to_cell_data()
dV = np.asarray(vc["Volume"]) * VOL_TO_SI                       # m^3
Um_tot = float(np.sum(np.asarray(vc["U_m"]) * dV))              # J  (time-averaged)
Ue_tot = float(np.sum(np.asarray(vc["U_e"]) * dV))              # J
imag_ratio = (np.median(np.linalg.norm(vol["B_imag"], axis=1))
              / np.median(np.linalg.norm(vol["B_real"], axis=1)))
print(f"U_m_tot={Um_tot:.5e} J  U_e_tot={Ue_tot:.5e} J  |B_imag|/|B_real|~{imag_ratio:.1e}")

# cross-check against Palace's own domain integral
try:
    dE = pd.read_csv(f"{OF}/domain-E.csv").to_numpy()
    emag = dE[CYC - 1, 2]
    print(f"cross-check: domain-E.csv E_mag = {emag:.5e} J   ratio U_m_tot/E_mag = {Um_tot/emag:.4f}")
except Exception as e:                                          # noqa
    emag = np.nan
    print(f"cross-check: domain-E.csv not read ({e})")

# horizontal cut Z_OFF_UM above the metal plane (mesh coords in mm)
sl = vol.slice(normal="z", origin=(0, 0, Z_OFF_UM * 1e-3)).compute_cell_sizes(area=True)
sc = sl.cell_centers().points                                   # mm
slc = sl.point_data_to_cell_data()
# eigenmode phase lives in the imaginary part here; keep B_z COMPLEX (Tesla)
sBz = np.asarray(slc["B_real"])[:, 2] + 1j * np.asarray(slc["B_imag"])[:, 2]
sA = np.asarray(sl["Area"]) * AREA_TO_SI                        # m^2

def bz_at(x_mm, y_mm, r_mm=0.004):
    m = (np.abs(sc[:, 0] - x_mm) < r_mm) & (np.abs(sc[:, 1] - y_mm) < r_mm)
    return complex(np.mean(sBz[m])) if m.any() else np.nan

x_g_mm = ((A_UM + B_UM) / 2) / 1e3
sqrt_fac = np.sqrt((((A_UM + B_UM) / 2) ** 2 - A_UM ** 2)
                   * (B_UM ** 2 - ((A_UM + B_UM) / 2) ** 2)) * 1e-12  # um^2 -> m^2

Phi = np.zeros(n_loops, complex)          # Wb  (mode-normalised, complex)
I_cpw = np.zeros(n_loops)                 # A   from the exact CPW gap formula
I_twin = np.zeros(n_loops)                # A   crude twin-lead
for k in range(n_loops):
    x0 = (ARRAY_X0_UM + k * PITCH_UM) / 1e3
    inx = (sc[:, 0] > x0 + 0.004) & (sc[:, 0] < x0 + 0.037) & (np.abs(sc[:, 1]) < 0.0125)
    Phi[k] = complex(np.sum(sBz[inx] * sA[inx]))
    xc = x0 + 0.020
    vals = []
    for leg in (-HAIRPIN_LEG_Y_UM, HAIRPIN_LEG_Y_UM):
        for side in (-1, 1):
            yg = (leg + side * (A_UM + B_UM) / 2) / 1e3
            bz = bz_at(xc, yg)
            if np.isfinite(bz):
                vals.append(abs(bz) * 4 * KK * sqrt_fac / (MU0 * (B_UM * 1e-6)))
    I_cpw[k] = np.mean(vals) if vals else np.nan
    bz0 = bz_at(xc, 0.0)
    I_twin[k] = abs(bz0) * np.pi * (HAIRPIN_LEG_Y_UM * 1e-6) / MU0

phi0g = np.angle(Phi[np.nanargmax(np.abs(Phi))])       # global mode phase
Phi_s = np.real(Phi * np.exp(-1j * phi0g))             # signed real flux amplitude (Wb)
Phi_mag = np.abs(Phi)
I_r = np.nanmax(I_cpw)                                 # modal current amplitude (peak, A)
I_r_twin = np.nanmax(I_twin)
L_r = 2.0 * Um_tot / I_r ** 2                          # lumped-equivalent inductance (H); Um_tot = peak E_mag
I_zpf = np.sqrt(HBAR * w_r / (2.0 * L_r))              # A
M_cpw = Phi_s / I_r                                    # H  (signed)
M_loc = Phi_mag / I_cpw                                # H  geometric per-loop mutual (expect ~flat)
c_k = np.pi * Phi_s * np.sqrt(HBAR * w_r / (4.0 * Um_tot)) / PHI0
# consistency: c_k should also equal pi * M_cpw * I_zpf / Phi0
c_k_viaM = np.pi * M_cpw * I_zpf / PHI0

print(f"\nL_r = {L_r*1e9:.4f} nH   I_r(peak) = {I_r:.4e} A   I_r^zpf = {I_zpf:.4e} A")
print(f"{'loop':>4} {'Phi_k (Wb)':>13} {'I_cpw (A)':>11} {'M_k (H)':>12} {'M_loc (H)':>12} "
      f"{'c_k':>11} {'c_k viaM':>11} {'p_EPR':>11}")
for k in range(n_loops):
    print(f"{k:>4} {Phi_s[k]:>13.3e} {I_cpw[k]:>11.3e} {M_cpw[k]:>12.3e} {M_loc[k]:>12.3e} "
          f"{c_k[k]:>11.3e} {c_k_viaM[k]:>11.3e} {ploop_csv[k]:>11.2e}")

# port-EPR consistency relation: p_loop,k = M_k^2 / (2 L_J L_r)
LJ_nH = float(OF.split("_Lj_")[1].split("nH")[0]) if "_Lj_" in OF else np.nan
if np.isfinite(LJ_nH):
    p_pred = M_cpw ** 2 / (2.0 * LJ_nH * 1e-9 * L_r)
    print(f"\np_loop predicted from M (p = M^2/(2 L_J L_r), L_J={LJ_nH:.0f} nH):")
    for k in (0, n_loops // 2, n_loops - 1):
        print(f"  loop {k:2d}:  p_pred = {p_pred[k]:.2e}   p_EPR = {ploop_csv[k]:.2e}   "
              f"ratio = {ploop_csv[k]/p_pred[k]:.1f}")
else:
    p_pred = np.full(n_loops, np.nan)

print(f"\nf_r = {f_r/1e9:.4f} GHz")
print(f"M_k     (H): {np.nanmin(M_cpw):.3e} .. {np.nanmax(M_cpw):.3e}   (Phi_s / I_r)")
print(f"M_loc   (H): mean {np.nanmean(M_loc):.3e} +/- {np.nanstd(M_loc):.2e}  "
      f"({100*np.nanstd(M_loc)/np.nanmean(M_loc):.0f}% spread -> flat = method clean)")
print(f"c_array (signed sum) = {np.nansum(c_k):.4e} rad")
print(f"c_array (sum |c_k|)  = {np.nansum(np.abs(c_k)):.4e} rad")
print(f"I_cpw peak / I_twin peak = {I_r/I_r_twin:.2f}  "
      f"(the two current models -- if this is far from 1 the current is uncertain)")

# ---- raw-data dump the user asked for ----
npz = os.path.join(os.path.dirname(OF.rstrip("/")),
                   f"M_bfield_cyc{CYC}_z{Z_OFF_UM:g}.npz")
np.savez(npz,
         loop=np.arange(n_loops),
         Phi_complex=Phi,                 # Wb, complex, mode-normalised
         Phi_signed=Phi_s,               # Wb, real, global-phase-aligned
         Phi_abs=Phi_mag,                # Wb
         I_cpw=I_cpw,                     # A, per loop, exact CPW gap formula
         I_twin=I_twin,                   # A, per loop, twin-lead model
         I_r_peak=I_r,                    # A, scalar, standing-wave antinode
         I_r_zpf=I_zpf,                   # A, scalar
         M_k=M_cpw,                        # H, per loop (Phi_s / I_r_peak)
         M_loc=M_loc,                      # H, per loop (|Phi| / I_local)
         c_k=c_k,                          # rad, per loop
         c_k_viaM=c_k_viaM,                # rad, per loop (cross-check route)
         p_EPR=ploop_csv,                  # Palace port-EPR per-loop participation
         p_pred_from_M=p_pred,             # M^2/(2 L_J L_r)
         L_r_H=L_r, U_m_tot_J=Um_tot, U_e_tot_J=Ue_tot,
         E_mag_domainE_J=emag,
         f_r_Hz=f_r, w_r=w_r, z_offset_um=Z_OFF_UM, cycle=CYC,
         L_J_nH=LJ_nH, n_loops=n_loops,
         A_um=A_UM, B_um=B_UM, Kk=KK, pitch_um=PITCH_UM, array_x0_um=ARRAY_X0_UM)
print(f"\nsaved raw arrays -> {npz}")

fig, ax = plt.subplots(4, 1, figsize=(11, 12), sharex=True)
ax[0].bar(range(n_loops), Phi_s * 1e15);        ax[0].set_title(f"Phi_k  (fWb)   z = {Z_OFF_UM:.1f} um above plane")
ax[1].bar(range(n_loops), I_cpw, color="#c85");  ax[1].set_title("local resonator current at each loop  (A, exact CPW gap formula)")
ax[2].bar(range(n_loops), M_loc * 1e12, color="#3a7"); ax[2].set_title("M_k^loc = Phi_k / I_local  (pH, geometric, expect ~flat)")
ax[3].bar(range(n_loops), np.abs(c_k), color="#47a"); ax[3].set_title("|c_k|  (rad)"); ax[3].set_xlabel("loop index")
for a in ax:
    a.grid(axis="y", alpha=.25)
plt.suptitle(f"{OF}\nc_array = {np.nansum(c_k):.3e} rad   M_k ~ {np.nanmean(M_cpw)*1e12:.2f} pH   "
             f"L_r = {L_r*1e9:.2f} nH   B_imag/B_real~{imag_ratio:.0e}")
plt.tight_layout()
plt.savefig("/tmp/M_bfield.png", dpi=115)
print("saved /tmp/M_bfield.png")
