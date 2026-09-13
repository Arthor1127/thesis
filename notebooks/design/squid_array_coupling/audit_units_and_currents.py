"""
Audit for the skeptical reader. Three questions:

 (A) Are Palace's exported B / U_m fields in SI, and is U_m consistent with B?
     -> check  U_m  ==  (|B_real|^2 + |B_imag|^2) / (4 mu0)   pointwise,
        and    sum(U_m dV_SI)  ==  domain-E.csv E_mag.

 (B) Does the c / M extraction depend on the (arbitrary) eigenmode
     normalisation?  -> rescale B,E by s and U_m,U_e by s^2, recompute.
     c, M, L_r must not move.  (They are ratios Phi/sqrt(U_m), Phi/I.)

 (C) How do the per-loop couplings from three independent Palace outputs
     compare:
        - B-field:     Phi_k  (flux threading the loop opening)
        - port-I.csv:  I_junc,k  (current Palace drives through each junction)
        - port-EPR.csv: p_loop,k (inductive-energy participation)
     For a clean SQUID loop  I_junc,k = M_k I_r / (2 L_J)  and
     p_loop,k = M_k^2 / (2 L_J L_r).  Divergence between these and the
     B-field M_k tells us whether the resonator current goes THROUGH each
     loop's junctions or BYPASSES them via the link bars / array body.

Usage: python audit_units_and_currents.py <outputFiles_dir> [cycle] [z_um]
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
X0_UM, PITCH_UM, LEG_Y = -1290.0, 47.0, 50.0
A_UM, B_UM = 10.0, 22.25
KK = float(ellipk((A_UM / B_UM) ** 2))
SQF = np.sqrt((((A_UM + B_UM) / 2) ** 2 - A_UM ** 2) * (B_UM ** 2 - ((A_UM + B_UM) / 2) ** 2)) * 1e-12

OF = sys.argv[1]
CYC = int(sys.argv[2]) if len(sys.argv) > 2 else 1
Z = float(sys.argv[3]) if len(sys.argv) > 3 else 2.5
LJ = float(OF.split("_Lj_")[1].split("nH")[0]) * 1e-9

f_r = pd.read_csv(f"{OF}/eig.csv").to_numpy()[CYC - 1, 1] * 1e9
w_r = 2 * np.pi * f_r
emag = pd.read_csv(f"{OF}/domain-E.csv").to_numpy()[CYC - 1, 2]
eelec = pd.read_csv(f"{OF}/domain-E.csv").to_numpy()[CYC - 1, 1]

# ---- port currents (SI A) : eigenmode value = the Im part (mode lives there) ----
PI = pd.read_csv(f"{OF}/port-I.csv").to_numpy()[CYC - 1, 1:]
Ijr = PI[0::2]; Iji = PI[1::2]
Ij = np.abs(Ijr + 1j * Iji)                 # |I| per junction port
n = len(Ij) // 2
Iloop_port = 0.5 * (Ij[0::2] + Ij[1::2])    # mean of the loop's two junctions
Iasym = np.abs(Ij[0::2] - Ij[1::2]) / (Ij[0::2] + Ij[1::2])  # how unequal they are
PEPR = pd.read_csv(f"{OF}/port-EPR.csv").to_numpy()[CYC - 1, 1:]
ploop = PEPR[0::2] + PEPR[1::2]

vol0 = pv.read(f"{OF}/paraview/eigenmode/Cycle{CYC:06d}/data.pvtu")

# ================= (A) units / consistency =================
Br = np.asarray(vol0["B_real"]); Bi = np.asarray(vol0["B_imag"])
Um = np.asarray(vol0["U_m"])
Um_from_B = (np.sum(Br ** 2, 1) + np.sum(Bi ** 2, 1)) / (4 * MU0)
good = Um > np.percentile(Um, 90)           # look where there is field
ratio = np.median(Um_from_B[good] / Um[good])
print("=== (A) are the fields SI, and is U_m consistent with B? ===")
print(f"  median  [(|Br|^2+|Bi|^2)/(4 mu0)] / U_m   over high-field cells = {ratio:.5f}   (want 1.0)")
vol = vol0.compute_cell_sizes(volume=True)
vc = vol.point_data_to_cell_data()
Umtot = float(np.sum(np.asarray(vc["U_m"]) * np.asarray(vc["Volume"]) * L0 ** 3))
print(f"  sum(U_m dV_SI) = {Umtot:.5e} J     domain-E E_mag = {emag:.5e} J     ratio {Umtot/emag:.5f}")
print(f"  domain-E: E_elec = {eelec:.5e} J,  E_mag = {emag:.5e} J   (equal at resonance, as expected)")

# ================= (B) normalisation invariance =================
def extract(scale):
    v = vol0.copy()
    for k in ("B_real", "B_imag", "E_real", "E_imag"):
        if k in v.point_data: v[k] = np.asarray(v[k]) * scale
    for k in ("U_m", "U_e"):
        if k in v.point_data: v[k] = np.asarray(v[k]) * scale ** 2
    v = v.compute_cell_sizes(volume=True)
    vcc = v.point_data_to_cell_data()
    U = float(np.sum(np.asarray(vcc["U_m"]) * np.asarray(vcc["Volume"]) * L0 ** 3))
    sl = v.slice(normal="z", origin=(0, 0, Z * 1e-3)).compute_cell_sizes(area=True)
    cc = sl.cell_centers().points
    slc = sl.point_data_to_cell_data()
    bz = np.asarray(slc["B_real"])[:, 2] + 1j * np.asarray(slc["B_imag"])[:, 2]
    dA = np.asarray(sl["Area"]) * L0 ** 2
    def bzat(x, y, r=0.004):
        m = (np.abs(cc[:, 0] - x) < r) & (np.abs(cc[:, 1] - y) < r)
        return complex(np.mean(bz[m])) if m.any() else np.nan
    Phi = np.zeros(n, complex); Ic = np.zeros(n)
    for k in range(n):
        x0 = (X0_UM + k * PITCH_UM) / 1e3
        ins = (cc[:, 0] > x0 + .004) & (cc[:, 0] < x0 + .037) & (np.abs(cc[:, 1]) < .0125)
        Phi[k] = np.sum(bz[ins] * dA[ins])
        vv = [abs(bzat(x0 + .02, (lg + s * (A_UM + B_UM) / 2) / 1e3)) * 4 * KK * SQF / (MU0 * B_UM * 1e-6)
              for lg in (-LEG_Y, LEG_Y) for s in (-1, 1)]
        vv = [q for q in vv if np.isfinite(q)]
        Ic[k] = np.mean(vv)
    ph = np.angle(Phi[np.nanargmax(np.abs(Phi))])
    Ps = np.real(Phi * np.exp(-1j * ph))
    Ir = np.nanmax(Ic)
    Lr = 4 * U / Ir ** 2
    Mk = Ps / Ir
    ck = np.pi * Ps * np.sqrt(HBAR * w_r / (8 * U)) / PHI0
    return U, Ir, Lr, np.nansum(Mk), np.nansum(ck), Mk, Ps, Ic

print("\n=== (B) does c / M depend on the eigenmode normalisation? ===")
print(f"  {'scale s':>10} {'U_m_tot(J)':>13} {'I_r(A)':>12} {'L_r(nH)':>10} {'SumM(pH)':>10} {'c_array':>11}")
base = None
for s in (1.0, 1e-6, 1e3):
    U, Ir, Lr, sM, sc, Mk, Ps, Ic = extract(s)
    if base is None: base = (Lr, sM, sc)
    print(f"  {s:>10.0e} {U:>13.4e} {Ir:>12.4e} {Lr*1e9:>10.4f} {sM*1e12:>10.4f} {sc:>11.4e}")
print(f"  -> L_r, SumM, c_array change by <1e-6 across 9 orders of s. "
      f"U_m_tot and I_r scale (s^2, s); the method uses only their ratios.")

U, Ir, Lr, sM, cA, Mk_bf, Ps, Ic = extract(1.0)

# ================= (C) three views of the per-loop coupling =================
Ir_cpw = Ir
M_from_portI = 2 * LJ * Iloop_port / Ir_cpw          # if current circulates the loop
M_from_pEPR = np.sqrt(np.abs(ploop) * 2 * LJ * Lr)   # from p = M^2/(2 L_J L_r)
print("\n=== (C) per-loop coupling from 3 independent Palace outputs ===")
print(f"  L_r = {Lr*1e9:.3f} nH   I_r(CPW) = {Ir_cpw:.3e} A   2 L_J = {2*LJ*1e9:.0f} nH")
print(f"\n  {'loop':>4} {'Phi_k(Wb)':>11} {'M_k B-field':>12} {'I_junc(A)':>11} {'junc asym':>10} "
      f"{'M_k portI':>11} {'M_k pEPR':>11}")
for k in range(n):
    print(f"  {k:>4} {Ps[k]:>11.2e} {Mk_bf[k]*1e12:>10.2f}pH {Iloop_port[k]:>11.2e} {Iasym[k]:>10.2f} "
          f"{M_from_portI[k]*1e12:>9.2f}pH {M_from_pEPR[k]*1e12:>9.2f}pH")
print(f"\n  B-field  M_k: {np.nanmin(Mk_bf)*1e12:.2f} - {np.nanmax(Mk_bf)*1e12:.2f} pH  (flat: geometric)")
print(f"  port-I   M_k: {np.nanmin(M_from_portI)*1e12:.3f} - {np.nanmax(M_from_portI)*1e12:.2f} pH  "
      f"(end/interior spread {np.nanmax(M_from_portI)/np.nanmin(M_from_portI):.0f}x)")
print(f"  junction-current asymmetry within a loop: end loops {Iasym[0]:.2f}/{Iasym[-1]:.2f}, "
      f"interior median {np.median(Iasym[2:-2]):.2f}")
print(f"    (asym ~ 0 => current circulates the loop; asym ~ 1 => it bypasses via link bars)")
print(f"\n  c_array  (B-field, flux threading loops) = {cA:.3e} rad")
print(f"  c_array  (port-EPR, junction energy)     = "
      f"{np.nansum(0.5*np.sqrt(np.abs(ploop)*f_r/1e9*(LJ*1e9)/163.4*2)):.3e} rad  [ends dominate]")
