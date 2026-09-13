"""
Cross-check the per-loop coupling from the eigenmode FIELDS against Palace's
port-EPR.csv. The array runs down a resonator hairpin; the question is
whether all loops couple (a monotonic ramp in U_m / J_s) or only the two
ends (what port-EPR.csv reports).

Usage: python analyze_loop_fields.py <run_outputFiles_dir> [cycle]
"""
import sys
import numpy as np
import pandas as pd
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pv.OFF_SCREEN = True
OF = sys.argv[1]
cyc = int(sys.argv[2]) if len(sys.argv) > 2 else 1

m = pv.read(f"{OF}/paraview/eigenmode_boundary/Cycle{cyc:06d}/data.pvtu")
m = m.compute_cell_sizes(area=True, length=False, volume=False)
mc = m.point_data_to_cell_data()
attr = np.asarray(m["attribute"])
area = np.asarray(m["Area"])
Um = np.asarray(mc["U_m"])
Ue = np.asarray(mc["U_e"])
Js = np.linalg.norm(mc["J_s_real"], axis=1)
cc = m.cell_centers().points

P = pd.read_csv(f"{OF}/port-EPR.csv").to_numpy()[:, 1:]
ploop_csv = P[cyc - 1][0::2] + P[cyc - 1][1::2]     # resonator mode = row cyc-1

Um_tot = float(np.sum(Um * area))
print(f"boundary Sum(U_m dA) = {Um_tot:.4e}   Sum(U_e dA) = {np.sum(Ue*area):.4e}")

# array metal is attribute 81 (verified from the mesh dump)
mx = cc[attr == 81, 0]
mJs = Js[attr == 81]
mArea = area[attr == 81]

rows = []
for k in range(23):
    xc = -1.290 + k * 0.047 + 0.0235
    sm = np.abs(mx - xc) < 0.023
    j_int = float(np.sum(mJs[sm] * mArea[sm]))          # circulating current proxy
    sp = (attr == 2 * k + 1) | (attr == 2 * k + 2)      # this loop's 2 port surfaces
    um_port = float(np.sum(Um[sp] * area[sp]) / Um_tot)  # participation from the field
    rows.append((k, um_port, j_int, ploop_csv[k]))

print(f"{'loop':>4} {'U_m_port/tot':>14} {'Js_on_metal':>14} {'port-EPR p':>14}")
for k, up, ji, pc in rows:
    print(f"{k:>4} {up:>14.3e} {ji:>14.3e} {pc:>14.3e}")

up = np.array([r[1] for r in rows])
ji = np.array([r[2] for r in rows])
pc = np.abs([r[3] for r in rows])

fig, ax = plt.subplots(3, 1, figsize=(11, 10))
ax[0].bar(range(23), up); ax[0].set_title("participation from field: integral of U_m over each loop's port surfaces / total")
ax[1].bar(range(23), ji, color="#a5522a"); ax[1].set_title("|J_s_real| on array metal near each loop  (circulating current)")
ax[2].bar(range(23), pc, color="#2266aa"); ax[2].set_yscale("log")
ax[2].set_title("|p_loop| from Palace port-EPR.csv")
for a in ax:
    a.set_xlabel("loop index"); a.set_xticks(range(0, 23, 2))
plt.suptitle(f"{OF}\nresonator mode, cycle {cyc}")
plt.tight_layout()
plt.savefig("/tmp/loop_compare.png", dpi=110)
print("saved /tmp/loop_compare.png")

# --- E across the junction gaps (what the lumped port actually senses) ---
Er = mc["E_real"]
print("\nloop  |E_real| integrated over the 2 junction-port surfaces:")
eL = []
for k in range(23):
    sp = (attr == 2 * k + 1) | (attr == 2 * k + 2)
    e = float(np.sum(np.linalg.norm(Er[sp], axis=1) * area[sp]))
    eL.append(e)
eL = np.array(eL)
for k in range(23):
    print(f"  {k:2d}  {eL[k]:.3e}   {'#'*int(50*eL[k]/eL.max())}")
