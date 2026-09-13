"""What HFSS eigenmode saw vs PALACE, for the 3 requested geometries.

Left : PALACE coupling Qc rises with the coupling gap, but the HFSS resonator
       mode's total Q sits flat on the dielectric limit (Qi ~ 92 000) for all
       three -- the R=50 ports add no measurable loss to it, so HFSS gives no
       coupling Qc.
Right: the ports DO work -- the ~16.5 GHz feedline mode is damped from the
       no-port dielectric limit down to Q ~ 15-18 k.
"""
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
QI = 92231.0

rows = []
for p in sorted(glob.glob(os.path.join(HERE, "eigen_qc_*_result.json"))):
    if p.endswith(".seeded"):
        continue
    d = json.load(open(p))
    pr = d["palace_reference"]
    pal = pr.get("eigenmode_Qc") or pr.get("eigenmode_Qc_2port_R50")
    rm = d["resonator_mode"]
    fl = [m for m in d["all_modes"]
          if m.get("f_ghz") and 15.5 < m["f_ghz"] < 17.5]
    flq = min((m["Q"] for m in fl), default=np.nan)
    rows.append(dict(geom=d["geometry"].replace("_", "\n"), palace=pal,
                     qtot=rm.get("Q_total_hfss"), fr=rm.get("fr_ghz"),
                     feedQ=flq))

rows.sort(key=lambda r: r["palace"] or 0)
x = np.arange(len(rows))
lab = [r["geom"] for r in rows]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))

ax[0].set_yscale("log")
ax[0].bar(x - 0.2, [r["palace"] for r in rows], 0.4, color="#c0392b",
          label="PALACE eigenmode Qc")
ax[0].bar(x + 0.2, [r["qtot"] for r in rows], 0.4, color="#1e8449",
          label="HFSS resonator Q_total")
ax[0].axhline(QI, ls="--", color="k", lw=1, label="dielectric Qi (tanδ=1.2e-5)")
ax[0].axhspan(82000, 85000, color="grey", alpha=.2,
              label="PALACE confirmed ceiling")
for i, r in enumerate(rows):
    ax[0].annotate(f"{r['palace']:,}", (i - 0.2, r["palace"]), ha="center",
                   va="bottom", fontsize=8)
    ax[0].annotate(f"{r['qtot']:,.0f}", (i + 0.2, r["qtot"]), ha="center",
                   va="bottom", fontsize=8)
ax[0].set_xticks(x); ax[0].set_xticklabels(lab, fontsize=8)
ax[0].set_ylabel("Q"); ax[0].set_ylim(5e3, 3e5)
ax[0].set_title("HFSS gives no coupling Qc:\nresonator Q_total pinned at the dielectric limit")
ax[0].legend(fontsize=7.5, loc="lower right")

ax[1].bar(x, [r["feedQ"] for r in rows], 0.5, color="#2471a3",
          label="HFSS ~16.5 GHz feedline mode, with R=50 ports")
ax[1].axhline(85000, ls="--", color="k", lw=1,
              label="same mode, no ports (dielectric limit)")
for i, r in enumerate(rows):
    ax[1].annotate(f"{r['feedQ']:,.0f}", (i, r["feedQ"]), ha="center",
                   va="bottom", fontsize=8)
ax[1].set_xticks(x); ax[1].set_xticklabels(lab, fontsize=8)
ax[1].set_ylabel("Q"); ax[1].set_yscale("log"); ax[1].set_ylim(5e3, 1.2e5)
ax[1].set_title("the ports DO work:\nfeedline modes damped ~85k → ~15-18k")
ax[1].legend(fontsize=7.5)

fig.tight_layout()
fig.savefig(os.path.join(HERE, "qc_vs_palace.png"), dpi=130)
print("wrote qc_vs_palace.png")
