"""
Summarise the single-circular-island capacitance sweep:
  - island_cap_summary.csv        : full data table
  - island_cap_summary.png        : C vs radius / gap / thickness, with the
                                    coax-edge formula fit + a physical model that
                                    actually describes the data, + the
                                    sim/formula ratio.

Run from scripts/python/island_cap/ after fetching out/results/*.json
(or point it at a directory):  python summarize_island.py [results_dir]
"""
import os
import sys
import glob
import json
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RESDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "out", "results")
EPS0 = 8.8541878128e-12

rows = []
for fn in sorted(glob.glob(f"{RESDIR}/*.json")):
    d = json.load(open(fn))
    d.pop("maxwell_matrix_F", None)
    for s in ("radius", "gap", "thickness"):
        if d["tag"].endswith(s):
            d["sweep"] = s
    rows.append(d)
if not rows:
    sys.exit(f"no result json in {RESDIR}")

eps_r = rows[0]["eps_r"]
eps_eff = EPS0 * (1.0 + eps_r) / 2.0        # coplanar half-space value

# ---- CSV ----
keys = ["tag", "sweep", "radius_um", "gap_um", "thickness_nm", "eps_r", "b_over_a",
        "C_island_ground_fF", "C_island_self_fF", "C_formula_coax_fF",
        "ratio_sim_over_formula", "dof", "mesh_elements", "solve_s"]
csv_path = os.path.join(HERE, "island_cap_summary.csv")
with open(csv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
    w.writeheader()
    for r in sorted(rows, key=lambda x: (x["sweep"], x["radius_um"], x["gap_um"], x["thickness_nm"])):
        w.writerow(r)
print("wrote", csv_path)


def coax(eps_e, a_um, g_um, t_nm):
    a = np.asarray(a_um) * 1e-6
    b = (np.asarray(a_um) + np.asarray(g_um)) * 1e-6
    t = np.asarray(t_nm) * 1e-9
    return 2.0 * np.pi * eps_e * t / np.log(b / a)


def grab(sw, xkey):
    pts = sorted([r for r in rows if r["sweep"] == sw], key=lambda r: r[xkey])
    x = np.array([p[xkey] for p in pts])
    a = np.array([p["radius_um"] for p in pts]) * 1e-6
    g = np.array([p["gap_um"] for p in pts]) * 1e-6
    t = np.array([p["thickness_nm"] for p in pts]) * 1e-9
    C = np.array([p["C_island_ground_fF"] for p in pts]) * 1e-15
    return pts, x, a, g, t, C


fig, ax = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)

# ---------- radius ----------
_, x, a, g, t, C = grab("radius", "radius_um")
xs = np.linspace(x.min(), x.max(), 200)
# coax (fit eps_eff)
bas = coax(1.0, x, g[0] * 1e6, t[0] * 1e9)
e_fit = np.sum(bas * C) / np.sum(bas * bas)
# physical model (NO free parameters): isolated disc 8*eps_eff*a with a
# logarithmic coplanar-gap enhancement,  C = 8*eps_eff*a*(1 + (1/pi)*ln(a/g)),
# eps_eff = eps0*(1+eps_r)/2.
def disc_model(a, g):
    return 8.0 * eps_eff * a * (1.0 + (1.0 / np.pi) * np.log(a / g))
ax[0, 0].plot(x, C * 1e15, "o", ms=7, label="Palace")
ax[0, 0].plot(xs, coax(e_fit, xs, g[0] * 1e6, t[0] * 1e9) * 1e15, "-",
              label=fr"coax fit, $\varepsilon_{{eff}}={e_fit/EPS0:.0f}\,\varepsilon_0$ (265$\times$ off)")
ax[0, 0].plot(xs, disc_model(xs * 1e-6, g[0]) * 1e15, "--", color="C2",
              label=r"$8\varepsilon_{eff}a[1+\frac{1}{\pi}\ln(a/g)]$  (no free params)")
ax[0, 0].plot(xs, 8 * eps_eff * xs * 1e-6 * 1e15, ":", color="gray", label=r"isolated disc $8\varepsilon_{eff}a$")
ax[0, 0].set(xlabel="island radius $a$ (µm)", ylabel="C (fF)", title="radius sweep  (g=10 µm, t=100 nm)")
ax[0, 0].grid(alpha=.3); ax[0, 0].legend(fontsize=8)

# ---------- gap ----------
_, x, a, g, t, C = grab("gap", "gap_um")
xs = np.linspace(x.min(), x.max(), 200)
bas = coax(1.0, a[0] * 1e6, x, t[0] * 1e9)
e_fit_g = np.sum(bas * C) / np.sum(bas * bas)
# C = A + B*ln(a/g)
A_ = np.vstack([np.ones_like(x), np.log(a[0] / (x * 1e-6))]).T
cg, *_ = np.linalg.lstsq(A_, C, rcond=None)
ax[0, 1].plot(x, C * 1e15, "o", ms=7, label="Palace")
ax[0, 1].plot(xs, coax(e_fit_g, a[0] * 1e6, xs, t[0] * 1e9) * 1e15, "-",
              label=fr"coax fit, $\varepsilon_{{eff}}={e_fit_g/EPS0:.0f}\,\varepsilon_0$")
ax[0, 1].plot(xs, disc_model(a[0], xs * 1e-6) * 1e15, "--", color="C2",
              label=r"$8\varepsilon_{eff}a[1+\frac{1}{\pi}\ln(a/g)]$  (no free params)")
ax[0, 1].plot(xs, (cg[0] + cg[1] * np.log(a[0] / (xs * 1e-6))) * 1e15, "-.", color="C4", lw=1,
              label=fr"$A + B\ln(a/g)$ fit  ($B={cg[1]*1e15:.0f}$ fF)")
ax[0, 1].axhline(8 * eps_eff * a[0] * 1e15, ls=":", color="gray", label=r"isolated disc $8\varepsilon_{eff}a$")
ax[0, 1].set(xlabel="gap $g$ (µm)", ylabel="C (fF)", title="gap sweep  (a=100 µm, t=100 nm)")
ax[0, 1].grid(alpha=.3); ax[0, 1].legend(fontsize=8)

# ---------- thickness ----------
_, x, a, g, t, C = grab("thickness", "thickness_nm")
xs = np.linspace(x.min(), x.max(), 200)
ax[1, 0].plot(x, C * 1e15, "o", ms=8, label="Palace")
ax[1, 0].axhline(C.mean() * 1e15, color="C2", ls="--", label=f"flat, C = {C.mean()*1e15:.1f} fF")
ax[1, 0].plot(xs, coax(eps_eff, a[0] * 1e6, g[0] * 1e6, xs) * 1e15, "-", color="C1",
              label=r"coax $\propto t$ (formula)")
ax[1, 0].set(xlabel="film thickness $t$ (nm)", ylabel="C (fF)",
             title="thickness sweep  (a=100 µm, g=10 µm)")
ax[1, 0].grid(alpha=.3); ax[1, 0].legend(fontsize=8)
ax[1, 0].text(0.5, 0.5, "quasi-3D: C is thickness-\nindependent (0.0% over 30-800 nm)",
              transform=ax[1, 0].transAxes, ha="center", fontsize=9,
              bbox=dict(boxstyle="round", fc="w", alpha=.8))

# ---------- sim / coax-formula ratio ----------
allr = sorted(rows, key=lambda r: (r["sweep"], r["radius_um"], r["gap_um"], r["thickness_nm"]))
lab = [r["tag"].split("_", 1)[1].rsplit("_", 1)[0] for r in allr]
rat = [r["ratio_sim_over_formula"] for r in allr]
col = {"radius": "C0", "gap": "C1", "thickness": "C2"}
ax[1, 1].barh(range(len(rat)), rat, color=[col[r["sweep"]] for r in allr])
ax[1, 1].set_yticks(range(len(rat))); ax[1, 1].set_yticklabels(lab, fontsize=7)
ax[1, 1].set(xlabel="C(Palace) / C(coax formula)",
             title="the coax-edge term is <0.5% of the real C")
ax[1, 1].axvline(1.0, color="k", lw=0.8)
ax[1, 1].grid(alpha=.3, axis="x")

fig.suptitle("Single circular island: island$\\to$ground capacitance\n"
             f"(Palace electrostatic, quasi-3D, $\\varepsilon_r$={eps_r}; order-1, ~15% mesh uncertainty)",
             fontsize=12)
out_png = os.path.join(HERE, "island_cap_summary.png")
fig.savefig(out_png, dpi=140)
print("wrote", out_png)

_rad = [r for r in rows if r["sweep"] == "radius"]
_gap = [r for r in rows if r["sweep"] == "gap"]
_dm = lambda r: disc_model(r["radius_um"] * 1e-6, r["gap_um"] * 1e-6) * 1e15
report = {
    "eps_r": eps_r,
    "eps_eff_coplanar_over_eps0": eps_eff / EPS0,
    "radius_sweep": {"coax_fit_eps_eff_over_eps0": e_fit / EPS0,
                     "disc_model_median_sim_over_model": float(np.median([r["C_island_ground_fF"] / _dm(r) for r in _rad]))},
    "gap_sweep": {"coax_fit_eps_eff_over_eps0": e_fit_g / EPS0,
                  "disc_model_median_sim_over_model": float(np.median([r["C_island_ground_fF"] / _dm(r) for r in _gap])),
                  "logfit_A_fF": cg[0] * 1e15, "logfit_B_fF": cg[1] * 1e15},
    "thickness_sweep": {"C_mean_fF": float(np.mean([r["C_island_ground_fF"] for r in rows if r["sweep"] == "thickness"])),
                        "C_spread_fF": float(np.ptp([r["C_island_ground_fF"] for r in rows if r["sweep"] == "thickness"]))},
    "median_sim_over_coax": float(np.median([r["ratio_sim_over_formula"] for r in rows])),
}
json.dump(report, open(os.path.join(HERE, "island_cap_summary.json"), "w"), indent=2)
print(json.dumps(report, indent=2))
