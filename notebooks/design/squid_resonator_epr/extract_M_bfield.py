"""
Direct B-field extraction of the resonator<->SQUID mutual inductance M and the
flux-coupling constant c, for the SINGLE-SQUID design (600 um arm), as an
independent cross-check of the port-EPR value in C_EXTRACTION_RESULTS.md
(c_inf = 2.23e-4 rad for l1=l2=600um).

Same method as ../squid_array_coupling/extract_M_bfield.py, reduced to one
loop.  See that file and ../squid_array_coupling/C_CALCULATION.md 2b for the
derivation and the units/normalisation notes.  Summary:

  Phi = integral of B_z over the SQUID loop interior, on a horizontal cut a
        few um ABOVE the metal plane (B.n = 0 on the plane by the PEC BC).
        B_z kept COMPLEX -- the eigenmode phase is in the imaginary part.

  I_r = resonator modal current amplitude, from its own CPW gap field via the
        exact zero-thickness-PEC result (Lahl & Woerdenweber APL 81 505 2002):
          B_z,gap(x) = mu0 I b / (4 K(k) sqrt((x^2-a^2)(b^2-x^2))),  a<|x|<b
        inverted at the gap midpoint.  a = cpw_width/2 = 10 um,
        b = a + cpw_gap = 22.25 um, k = a/b, K(k) = 1.66.

  M   = Phi / I_r
  L_r = 2 U_m_tot / I_r^2         (U_m_tot = sum(U_m * dV_SI) = domain-E E_mag,
                                   the PEAK inductive energy 1/2 L_r I_r^2)
  c   = pi M I_r^zpf / Phi0 = pi Phi sqrt(hbar w_r / (4 U_m_tot)) / Phi0
        (I_r cancels; c needs only Phi and U_m_tot, both ~ eigenmode scale,
         ratio is scale-invariant)

Geometry (from build_design.py + circle_transmon_squid.py source, l1=l2=600um):
  - resonator hairpin fold: legs at y = +/- 50 um, x from -1330 to -150 um
    (same as the squid array)
  - Each SQUID branch: pad chord -> feed of length l1 -> JUNCTION -> tail of
    length l2 -> INTO GROUND (the tails do NOT rejoin; each goes to ground).
  - The SQUID LOOP for the fluxoid (loop_opening in the source, "the actual
    SQUID loop area") is bounded by junction a, the l1 feed of arm a, the pad,
    and the l1 feed of arm b: x in [-689 (junction), -89 (pad)], |y| < 8.7 um.
    The region BEYOND the junctions (upper_void, x < -689) is a ground-plane
    void, NOT part of the electrical loop -- do NOT integrate it.
  - loop interior to integrate: x in [-686, -92] um, |y| < 8.5 um.

port-EPR.csv columns for this design: p[1] = transmon, p[2] = SQUID probe arm.

Usage:  python extract_M_bfield.py <outputFiles_dir> [cycle] [z_offset_um]
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
AREA_TO_SI = L0 ** 2
VOL_TO_SI = L0 ** 3

CPW_W_UM = 20.0
CPW_GAP_UM = 12.25
HAIRPIN_LEG_Y_UM = 50.0
A_UM = CPW_W_UM / 2               # 10
B_UM = A_UM + CPW_GAP_UM          # 22.25
KMOD = A_UM / B_UM
KK = float(ellipk(KMOD ** 2))

# SQUID loop = loop_opening: pad chord (x ~ -89) to the junctions (x ~ -689).
LOOP_X_LO, LOOP_X_HI = -686.0, -92.0
LOOP_Y_HALF = 8.5
JUNCTION_X_UM = -690.0
# for reference: the full arm span and the beyond-junction void
FULL_X_LO, VOID_X_LO, VOID_X_HI = -1275.0, -1286.0, -692.0

OF = sys.argv[1]
CYC = int(sys.argv[2]) if len(sys.argv) > 2 else 1
Z_OFF_UM = float(sys.argv[3]) if len(sys.argv) > 3 else 2.5

f_r = pd.read_csv(f"{OF}/eig.csv").to_numpy()[CYC - 1, 1] * 1e9
w_r = 2 * np.pi * f_r
P = pd.read_csv(f"{OF}/port-EPR.csv").to_numpy()[:, 1:]
p_transmon = P[CYC - 1, 0]
p_probe = P[CYC - 1, 1]
try:
    L_probe_nH = float(OF.split("_Lprobe_")[1].split("_")[0])
except Exception:
    L_probe_nH = np.nan
print(f"cycle {CYC}: f_r = {f_r/1e9:.4f} GHz   p_transmon={p_transmon:.3e}  p_probe={p_probe:.3e}"
      f"   L_probe={L_probe_nH} nH   K(k={KMOD:.3f})={KK:.4f}")

# reference port-EPR c
if np.isfinite(L_probe_nH):
    E_probe_h = 163.4 / L_probe_nH                         # GHz
    c_portEPR = 0.5 * np.sqrt(abs(p_probe) * (f_r / 1e9) / (2 * E_probe_h))
    print(f"port-EPR:  c = 0.5 sqrt(p_probe f / (2 E_probe/h)) = {c_portEPR:.4e} rad "
          f"(sign of p_probe: {np.sign(p_probe):+.0f})")
else:
    c_portEPR = np.nan

vol = pv.read(f"{OF}/paraview/eigenmode/Cycle{CYC:06d}/data.pvtu").compute_cell_sizes(volume=True)
vc = vol.point_data_to_cell_data()
dV = np.asarray(vc["Volume"]) * VOL_TO_SI
Um_tot = float(np.sum(np.asarray(vc["U_m"]) * dV))
Ue_tot = float(np.sum(np.asarray(vc["U_e"]) * dV))
imag_ratio = (np.median(np.linalg.norm(vol["B_imag"], axis=1))
              / np.median(np.linalg.norm(vol["B_real"], axis=1)))
try:
    emag = pd.read_csv(f"{OF}/domain-E.csv").to_numpy()[CYC - 1, 2]
except Exception:
    emag = np.nan
print(f"U_m_tot={Um_tot:.5e} J   U_e_tot={Ue_tot:.5e} J   |B_imag|/|B_real|~{imag_ratio:.1e}")
print(f"cross-check: domain-E E_mag = {emag:.5e} J   ratio U_m_tot/E_mag = {Um_tot/emag:.5f}")

# units consistency: Palace U_m density should be (|Br|^2+|Bi|^2)/(2 mu0)
Br = np.asarray(vol["B_real"]); Bi = np.asarray(vol["B_imag"]); Um = np.asarray(vol["U_m"])
g = Um > np.percentile(Um, 90)
print(f"units: median [(|Br|^2+|Bi|^2)/(4 mu0)] / U_m over high-field cells = "
      f"{np.median((np.sum(Br[g]**2,1)+np.sum(Bi[g]**2,1))/(4*MU0)/Um[g]):.4f}  (0.5 => Palace U_m is peak)")

def _slice(z_um):
    s = vol.slice(normal="z", origin=(0, 0, z_um * 1e-3)).compute_cell_sizes(area=True)
    cc = s.cell_centers().points
    sc_ = s.point_data_to_cell_data()
    bz = np.asarray(sc_["B_real"])[:, 2] + 1j * np.asarray(sc_["B_imag"])[:, 2]
    dA = np.asarray(s["Area"]) * AREA_TO_SI
    return cc, bz, dA

def flux_region(z_um, x_lo, x_hi):
    cc, bz, dA = _slice(z_um)
    m = ((cc[:, 0] * 1e3 > x_lo) & (cc[:, 0] * 1e3 < x_hi)
         & (np.abs(cc[:, 1] * 1e3) < LOOP_Y_HALF))
    return complex(np.sum(bz[m] * dA[m]))

def flux_at_z(z_um):
    return flux_region(z_um, LOOP_X_LO, LOOP_X_HI)

# how the flux splits: SQUID loop (loop_opening) vs the beyond-junction void
for z in (0.3, 1.0):
    fo = abs(flux_region(z, LOOP_X_LO, LOOP_X_HI))
    fv = abs(flux_region(z, VOID_X_LO, VOID_X_HI))
    ff = abs(flux_region(z, FULL_X_LO, LOOP_X_HI))
    print(f"z={z} um:  |Phi| loop_opening={fo:.3e}   beyond-junction void={fv:.3e}   "
          f"full arm span={ff:.3e}")

# The loop is a THIN slot (~17 um), so B_z there rises with height as leg flux
# curves through. The fluxoid Phi = closed integral of A.dl around the z=0
# superconducting path = integral of B_z over the slot at z->0+.  Sample low
# and extrapolate.
z_scan = np.array([0.5, 1.0, 1.5, 1.75, 2.0, 2.25, 2.5, 3.0, 4.0, 5.0])
phi_scan = np.array([abs(flux_at_z(z)) for z in z_scan])
# The near-surface points under-read (cells straddle the zero-thickness PEC
# sheets); the far points over-read (they catch resonator-leg flux curving
# through above the loop). The plateau in between is the loop-linking flux.
mid = (z_scan >= 1.5) & (z_scan <= 3.0)
phi0_lin = float(np.mean(phi_scan[mid]))
phi0_quad = float(np.median(phi_scan[mid]))
print("\nPhi(z) over the SQUID loop (|Phi|, Wb):")
for z, p in zip(z_scan, phi_scan):
    print(f"  z={z:4.2f} um : {p:.4e}")
print(f"  plateau (z=1.5-3 um):  mean {phi0_lin:.4e}   median {phi0_quad:.4e}")

# horizontal cut at the requested Z_OFF_UM (for the maps / detailed checks)
sl = vol.slice(normal="z", origin=(0, 0, Z_OFF_UM * 1e-3)).compute_cell_sizes(area=True)
sc = sl.cell_centers().points
slc = sl.point_data_to_cell_data()
sBz = np.asarray(slc["B_real"])[:, 2] + 1j * np.asarray(slc["B_imag"])[:, 2]
sA = np.asarray(sl["Area"]) * AREA_TO_SI

def bz_at(x_mm, y_mm, r_mm=0.004):
    m = (np.abs(sc[:, 0] - x_mm) < r_mm) & (np.abs(sc[:, 1] - y_mm) < r_mm)
    return complex(np.mean(sBz[m])) if m.any() else np.nan

# ---- flux through the loop interior ----
inx = ((sc[:, 0] * 1e3 > LOOP_X_LO) & (sc[:, 0] * 1e3 < LOOP_X_HI)
       & (np.abs(sc[:, 1] * 1e3) < LOOP_Y_HALF))
Phi = complex(np.sum(sBz[inx] * sA[inx]))
phi_ph = np.angle(Phi)
Phi_s = float(np.real(Phi * np.exp(-1j * phi_ph)))          # = |Phi|, phase-aligned
Phi_mag = abs(Phi)
print(f"\nloop interior: {inx.sum()} cells, area = {sA[inx].sum()*1e12:.1f} um^2")
print(f"Phi = {Phi:.4e} Wb   |Phi| = {Phi_mag:.4e}   Phi_signed = {Phi_s:.4e}")

# robustness: exclude a +/- band around the junction x
for band_um in (20, 40, 80):
    m2 = inx & (np.abs(sc[:, 0] * 1e3 - JUNCTION_X_UM) > band_um)
    print(f"  drop +/-{band_um:>2} um around junction (x={JUNCTION_X_UM:.0f}): "
          f"|Phi| = {abs(np.sum(sBz[m2]*sA[m2])):.4e}  ({abs(np.sum(sBz[m2]*sA[m2]))/Phi_mag*100:.1f}%)")
# robustness: drop the brightest 1% of cells
thr = np.percentile(np.abs(sBz[inx]), 99)
m3 = inx & (np.abs(sBz) <= thr)
print(f"  drop brightest 1% of cells: |Phi| = {abs(np.sum(sBz[m3]*sA[m3])):.4e}  "
      f"({abs(np.sum(sBz[m3]*sA[m3]))/Phi_mag*100:.1f}%)   max/median |Bz| in window = "
      f"{np.abs(sBz[inx]).max()/np.median(np.abs(sBz[inx])):.1f}")

# ---- resonator current along the fold, from the CPW gap field ----
sqrt_fac = np.sqrt((((A_UM + B_UM) / 2) ** 2 - A_UM ** 2)
                   * (B_UM ** 2 - ((A_UM + B_UM) / 2) ** 2)) * 1e-12   # um^2 -> m^2
xs_um = np.arange(LOOP_X_LO + 30, LOOP_X_HI - 30, 40.0)
I_cpw = np.full(xs_um.size, np.nan)
for i, xu in enumerate(xs_um):
    vals = []
    for leg in (-HAIRPIN_LEG_Y_UM, HAIRPIN_LEG_Y_UM):
        for side in (-1, 1):
            yg = (leg + side * (A_UM + B_UM) / 2) / 1e3
            bz = bz_at(xu / 1e3, yg)
            if np.isfinite(bz):
                vals.append(abs(bz) * 4 * KK * sqrt_fac / (MU0 * (B_UM * 1e-6)))
    if vals:
        I_cpw[i] = np.mean(vals)
I_r = np.nanmax(I_cpw)                                       # antinode amplitude
I_at_loop = np.nanmean(I_cpw)                                # avg over the loop span

L_r = 2.0 * Um_tot / I_r ** 2
I_zpf = np.sqrt(HBAR * w_r / (2.0 * L_r))

def c_and_M(phi):
    M_ = phi / I_r
    c_ = np.pi * phi * np.sqrt(HBAR * w_r / (4.0 * Um_tot)) / PHI0
    return c_, M_

c_z, M_z = c_and_M(Phi_s)                       # at the requested Z_OFF_UM
c_lin, M_lin = c_and_M(phi0_lin)               # extrapolated to z=0 (linear)
c_quad, M_quad = c_and_M(phi0_quad)            # extrapolated to z=0 (quadratic)

print(f"\nI_cpw along the loop: min {np.nanmin(I_cpw):.3e}  max(antinode) {I_r:.3e}  mean {I_at_loop:.3e} A")
print(f"L_r = {L_r*1e9:.4f} nH   I_r^zpf = {I_zpf:.4e} A")
print(f"\n{'flux surface':>22} {'Phi (Wb)':>12} {'M (pH)':>9} {'c (rad)':>11} {'c/c_portEPR':>12}")
for lbl, phi, cc, mm in (("z = %.1f um" % Z_OFF_UM, Phi_s, c_z, M_z),
                         ("z -> 0 (linear)", phi0_lin, c_lin, M_lin),
                         ("z -> 0 (quadratic)", phi0_quad, c_quad, M_quad)):
    r = cc / c_portEPR if np.isfinite(c_portEPR) else np.nan
    print(f"{lbl:>22} {phi:>12.3e} {mm*1e12:>9.2f} {cc:>11.3e} {r:>12.2f}")
print(f"{'port-EPR reference':>22} {'':>12} {'':>9} {c_portEPR:>11.3e} {1.0:>12.2f}")

if np.isfinite(L_probe_nH):
    p_pred = (phi0_quad / I_r) ** 2 / (L_probe_nH * 1e-9 * L_r)
    print(f"\np_probe predicted from M(z->0) (= M^2/(L_probe L_r)) = {p_pred:.3e}   "
          f"vs port-EPR p_probe = {p_probe:.3e}   ratio = {p_probe/p_pred:.2f}")
c_bf = c_quad; M = M_quad; c_bf_viaM = c_quad; M_local = M_lin

npz = os.path.join(os.path.dirname(OF.rstrip("/")), f"M_bfield_cyc{CYC}_z{Z_OFF_UM:g}.npz")
np.savez(npz, Phi=Phi, Phi_signed=Phi_s, I_cpw=I_cpw, xs_um=xs_um, I_r_peak=I_r,
         I_at_loop=I_at_loop, I_r_zpf=I_zpf, M_H=M, M_local_H=M_local, L_r_H=L_r,
         c_bfield=c_bf, c_bfield_viaM=c_bf_viaM, c_portEPR=c_portEPR,
         p_probe=p_probe, p_transmon=p_transmon, U_m_tot_J=Um_tot, E_mag_domainE_J=emag,
         f_r_Hz=f_r, L_probe_nH=L_probe_nH, z_offset_um=Z_OFF_UM,
         loop_x=(LOOP_X_LO, LOOP_X_HI), loop_y_half=LOOP_Y_HALF)
print(f"\nsaved -> {npz}")

fig, ax = plt.subplots(1, 2, figsize=(13, 4))
ax[0].plot(xs_um, I_cpw * 1e3, "o-", color="#c85")
ax[0].set_xlabel("x along the fold (um)"); ax[0].set_ylabel("I_cpw (mA)")
ax[0].set_title(f"resonator current over the SQUID span   f_r={f_r/1e9:.3f} GHz")
ax[0].grid(alpha=.3)
xw, yw, bw = sc[inx, 0] * 1e3, sc[inx, 1] * 1e3, np.real(sBz[inx] * np.exp(-1j * phi_ph))
sccat = ax[1].scatter(xw, yw, c=bw, s=6, cmap="RdBu_r",
                      vmin=-np.percentile(np.abs(bw), 98), vmax=np.percentile(np.abs(bw), 98))
ax[1].set_title("B_z in the loop interior (phase-aligned)")
ax[1].set_xlabel("x (um)"); ax[1].set_ylabel("y (um)")
plt.colorbar(sccat, ax=ax[1])
plt.suptitle(f"{OF}\nc_bfield = {c_bf:.3e}   c_portEPR = {c_portEPR:.3e}   M = {M*1e12:.2f} pH   L_r = {L_r*1e9:.2f} nH")
plt.tight_layout()
plt.savefig("/tmp/M_bfield_single.png", dpi=115)
print("saved /tmp/M_bfield_single.png")
