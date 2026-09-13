"""Qc vs R -- analytic, PALACE eigenmode, HFSS eigenmode."""
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

R, Qa, Qp, Qh, Qd = [], [], [], [], []
for p in sorted(glob.glob(os.path.join(HERE, "bareres_eig_R*_result.json"))):
    d = json.load(open(p))
    tag = d["tag"].replace("_eig", "")
    R.append(d["R_ohm"]); Qa.append(d["analytic_Qc"])
    Qp.append(d["palace_eigenmode_Qc"]); Qh.append(d["Qc_eigenmode"])
    dp = os.path.join(HERE, tag + "_result.json")
    q = None
    if os.path.exists(dp):
        dd = json.load(open(dp))
        q = dd.get("Qc_oneport_fit")
        # sanity: driven fit must be within 3x of analytic to count
        if q and not (0.3 < q / d["analytic_Qc"] < 3):
            q = None
    Qd.append(q)
R = np.array(R)
Qd_ok = np.array([q if q else np.nan for q in Qd], float)

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].loglog(R, Qa, "k--o", label="analytic  (π/4)(R/Z0)")
ax[0].loglog(R, Qp, "s-", color="#c0392b", label="PALACE eigenmode")
ax[0].loglog(R, Qd_ok, "^-", color="#1e8449", label="HFSS driven S11 fit (this work)")
ax[0].loglog(R, Qh, "D:", color="#2471a3", alpha=.6, label="HFSS eigenmode (secondary)")
ax[0].set_xlabel("R (Ω)"); ax[0].set_ylabel("Qc")
ax[0].set_title("Qc vs R"); ax[0].grid(alpha=.3, which="both"); ax[0].legend(fontsize=8)

ax[1].semilogx(R, np.array(Qp) / np.array(Qa), "s-", color="#c0392b", label="PALACE eig / analytic")
ax[1].semilogx(R, Qd_ok / np.array(Qa), "^-", color="#1e8449", label="HFSS driven / analytic")
ax[1].semilogx(R, np.array(Qh) / np.array(Qa), "D:", color="#2471a3", alpha=.6, label="HFSS eig / analytic")
ax[1].axhline(1.0, color="k", ls="--", lw=1)
ax[1].set_xlabel("R (Ω)"); ax[1].set_ylabel("Qc / analytic")
ax[1].set_title("ratio to analytic — driven S11 flat at 1.20; eigenmodes drift")
ax[1].grid(alpha=.3); ax[1].legend(fontsize=8); ax[1].set_ylim(0.5, 1.4)

fig.tight_layout()
fig.savefig(os.path.join(HERE, "qc_vs_R.png"), dpi=130)
print("wrote qc_vs_R.png")
