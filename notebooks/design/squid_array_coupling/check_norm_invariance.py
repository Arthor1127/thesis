"""
Answer the question: does the c / M extraction depend on the ABSOLUTE scale of
Palace's eigenmode fields?

An eigenmode is defined only up to an overall complex scalar. Palace fixes one
(domain-E.csv shows every mode normalised to the same total energy ~1.18e-11 J),
but that constant is arbitrary. This script re-scales B, E, U_m, U_e by an
arbitrary factor and re-runs the extraction, to show c, M, L_r are invariant.

  B -> s B ,  E -> s E ,  U_m -> s^2 U_m ,  U_e -> s^2 U_e

  Phi_k   ~ s          (integral of B)
  I_cpw   ~ s          (linear in B, via the CPW gap formula)
  U_m_tot ~ s^2
  L_r = 4 U_m_tot / I_r^2                 ~ s^2 / s^2 = s^0   INVARIANT
  M_k = Phi_k / I_r                       ~ s   / s   = s^0   INVARIANT
  c_k = pi Phi_k sqrt(hbar w / 8 U_m_tot)/Phi0  ~ s / s = s^0 INVARIANT

Usage: python check_norm_invariance.py <outputFiles_dir> [cycle] [z_um]
"""
import sys
import numpy as np
import pandas as pd
import pyvista as pv
from scipy.special import ellipk

pv.OFF_SCREEN = True
MU0 = 4e-7 * np.pi
HBAR = 1.054571817e-34
PHI0 = 2.067833848e-15
L0 = 1e-3
ARRAY_X0_UM, PITCH_UM, LEG_Y_UM = -1290.0, 47.0, 50.0
A_UM, B_UM = 10.0, 22.25
KK = float(ellipk((A_UM / B_UM) ** 2))

OF = sys.argv[1]
CYC = int(sys.argv[2]) if len(sys.argv) > 2 else 1
Z = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0

f_r = pd.read_csv(f"{OF}/eig.csv").to_numpy()[CYC - 1, 1] * 1e9
w_r = 2 * np.pi * f_r
vol0 = pv.read(f"{OF}/paraview/eigenmode/Cycle{CYC:06d}/data.pvtu")
emag = pd.read_csv(f"{OF}/domain-E.csv").to_numpy()[CYC - 1, 2]
sqrt_fac = np.sqrt((((A_UM + B_UM) / 2) ** 2 - A_UM ** 2)
                   * (B_UM ** 2 - ((A_UM + B_UM) / 2) ** 2)) * 1e-12


def extract(scale):
    vol = vol0.copy()
    for k in ("B_real", "B_imag", "E_real", "E_imag"):
        if k in vol.point_data:
            vol[k] = np.asarray(vol[k]) * scale
    for k in ("U_m", "U_e"):
        if k in vol.point_data:
            vol[k] = np.asarray(vol[k]) * scale ** 2
    vol = vol.compute_cell_sizes(volume=True)
    vc = vol.point_data_to_cell_data()
    Um_tot = float(np.sum(np.asarray(vc["U_m"]) * np.asarray(vc["Volume"]) * L0 ** 3))

    sl = vol.slice(normal="z", origin=(0, 0, Z * 1e-3)).compute_cell_sizes(area=True)
    scc = sl.cell_centers().points
    slc = sl.point_data_to_cell_data()
    sBz = np.asarray(slc["B_real"])[:, 2] + 1j * np.asarray(slc["B_imag"])[:, 2]
    sA = np.asarray(sl["Area"]) * L0 ** 2

    def bz_at(x, y, r=0.004):
        m = (np.abs(scc[:, 0] - x) < r) & (np.abs(scc[:, 1] - y) < r)
        return complex(np.mean(sBz[m])) if m.any() else np.nan

    n = 23 if "nsq_23" in OF else int(OF.split("nsq_")[1].split("_")[0])
    Phi = np.zeros(n, complex)
    Icpw = np.zeros(n)
    for k in range(n):
        x0 = (ARRAY_X0_UM + k * PITCH_UM) / 1e3
        inx = (scc[:, 0] > x0 + 0.004) & (scc[:, 0] < x0 + 0.037) & (np.abs(scc[:, 1]) < 0.0125)
        Phi[k] = complex(np.sum(sBz[inx] * sA[inx]))
        vals = []
        for leg in (-LEG_Y_UM, LEG_Y_UM):
            for s in (-1, 1):
                yg = (leg + s * (A_UM + B_UM) / 2) / 1e3
                bz = bz_at(x0 + 0.020, yg)
                if np.isfinite(bz):
                    vals.append(abs(bz) * 4 * KK * sqrt_fac / (MU0 * (B_UM * 1e-6)))
        Icpw[k] = np.mean(vals)
    phi0g = np.angle(Phi[np.nanargmax(np.abs(Phi))])
    Phi_s = np.real(Phi * np.exp(-1j * phi0g))
    I_r = np.nanmax(Icpw)
    L_r = 4 * Um_tot / I_r ** 2
    M_k = Phi_s / I_r
    c_k = np.pi * Phi_s * np.sqrt(HBAR * w_r / (8 * Um_tot)) / PHI0
    return Um_tot, I_r, L_r, np.nansum(M_k), np.nansum(c_k), M_k[0], c_k[0]


print(f"{OF}\nf_r = {f_r/1e9:.4f} GHz   domain-E E_mag = {emag:.4e} J\n")
print(f"{'scale':>10} {'U_m_tot(J)':>13} {'I_r(a.u.)':>12} {'L_r(nH)':>10} "
      f"{'SumM(pH)':>10} {'c_array':>11} {'M_0(pH)':>10} {'c_0':>11}")
for s in (1.0, 1e-3, 1e3, 1e6, 3.7):
    Um, Ir, Lr, sM, sc, M0, c0 = extract(s)
    print(f"{s:>10.1e} {Um:>13.4e} {Ir:>12.4e} {Lr*1e9:>10.4f} "
          f"{sM*1e12:>10.4f} {sc:>11.4e} {M0*1e12:>10.4f} {c0:>11.4e}")

print("\n-> U_m_tot and I_r scale with the eigenmode normalisation;")
print("   L_r, M, c do NOT.  The method uses only ratios (Phi/sqrt(U_m), Phi/I_r).")
print("   What it DOES need right: the mesh length unit (L0=1e-3) converted to SI,")
print("   and U_m consistent with B (Palace: U_m = |B|^2/4mu0).  Check vs domain-E:")
Um1 = extract(1.0)[0]
print(f"   sum(U_m dV_SI) = {Um1:.5e} J   vs domain-E E_mag = {emag:.5e} J   ratio {Um1/emag:.4f}")
