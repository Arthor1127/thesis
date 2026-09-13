"""
Collate all out/results/*.json island-capacitance runs into a CSV and produce
one fit figure per swept parameter, comparing the simulated island->ground
capacitance against the coaxial-edge model

    C = 2*pi * eps_eff * t / ln(b/a),   b = a + g

Fit: eps_eff is the single free parameter (a multiplicative constant); the
result is compared to the coplanar half-space value eps0*(1+eps_r)/2.

Usage:  python aggregate_island.py [results_dir]   (default: ./out/results)
"""
import os
import sys
import json
import glob

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESDIR = sys.argv[1] if len(sys.argv) > 1 else "./out/results"
OUT = os.path.dirname(RESDIR.rstrip("/")) or "."
EPS0 = 8.8541878128e-12

rows = []
for fn in sorted(glob.glob(f"{RESDIR}/*.json")):
    with open(fn) as f:
        d = json.load(f)
    d.pop("maxwell_matrix_F", None)
    # sweep label is embedded in the tag ("..._radius"/"_gap"/"_thickness")
    for s in ("radius", "gap", "thickness"):
        if d["tag"].endswith(s):
            d["sweep"] = s
    rows.append(d)

if not rows:
    sys.exit(f"no result json files in {RESDIR}")

import csv
keys = ["tag", "sweep", "radius_um", "gap_um", "thickness_nm", "eps_r", "b_over_a",
        "C_island_ground_fF", "C_island_self_fF", "C_formula_coax_fF",
        "ratio_sim_over_formula", "dof", "mesh_elements", "solve_s", "wall_s"]
csv_path = f"{OUT}/island_cap_summary.csv"
with open(csv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
    w.writeheader()
    for r in sorted(rows, key=lambda x: (x.get("sweep", ""), x["radius_um"], x["gap_um"], x["thickness_nm"])):
        w.writerow(r)
print(f"wrote {csv_path}  ({len(rows)} runs)")


def coax_C(eps_eff, a_um, g_um, t_nm):
    a = np.asarray(a_um) * 1e-6
    b = (np.asarray(a_um) + np.asarray(g_um)) * 1e-6
    t = np.asarray(t_nm) * 1e-9
    return 2.0 * np.pi * eps_eff * t / np.log(b / a)


eps_r = rows[0]["eps_r"]
eps_eff_coplanar = EPS0 * (1.0 + eps_r) / 2.0

sweeps = {
    "radius":    ("radius_um",    "island radius $a$ (um)"),
    "gap":       ("gap_um",       "gap $g$ (um)"),
    "thickness": ("thickness_nm", "film thickness $t$ (nm)"),
}

fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
fit_report = {}

for ax, (sw, (xkey, xlabel)) in zip(axes, sweeps.items()):
    pts = sorted([r for r in rows if r.get("sweep") == sw], key=lambda r: r[xkey])
    if not pts:
        ax.set_title(f"{sw}: no data"); continue
    x = np.array([p[xkey] for p in pts])
    a = np.array([p["radius_um"] for p in pts])
    g = np.array([p["gap_um"] for p in pts])
    t = np.array([p["thickness_nm"] for p in pts])
    C = np.array([p["C_island_ground_fF"] for p in pts]) * 1e-15   # F

    # single-parameter fit: C = eps_eff * [ 2*pi*t/ln(b/a) ]
    basis = coax_C(1.0, a, g, t)            # C per unit eps_eff
    eps_fit = float(np.sum(basis * C) / np.sum(basis * basis))
    C_fit = eps_fit * basis
    C_nofit = coax_C(eps_eff_coplanar, a, g, t)
    ss_res = np.sum((C - C_fit) ** 2)
    ss_tot = np.sum((C - C.mean()) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    fit_report[sw] = {
        "eps_eff_fit": eps_fit,
        "eps_eff_fit_over_eps0": eps_fit / EPS0,
        "eps_eff_coplanar_over_eps0": eps_eff_coplanar / EPS0,
        "fit_over_coplanar": eps_fit / eps_eff_coplanar,
        "R2": r2,
        "median_sim_over_formula": float(np.median(C / C_nofit)),
    }

    xs = np.linspace(x.min(), x.max(), 200)
    # dense curve: hold the other two params at their (constant) sweep value
    a_d = np.full_like(xs, a[0]) if sw != "radius" else xs
    g_d = np.full_like(xs, g[0]) if sw != "gap" else xs
    t_d = np.full_like(xs, t[0]) if sw != "thickness" else xs

    ax.plot(x, C * 1e15, "o", ms=6, label="Palace (island$\\to$ground)")
    ax.plot(xs, coax_C(eps_fit, a_d, g_d, t_d) * 1e15, "-",
            label=f"coax fit  $\\varepsilon_{{eff}}/\\varepsilon_0$={eps_fit/EPS0:.2f}  ($R^2$={r2:.3f})")
    ax.plot(xs, coax_C(eps_eff_coplanar, a_d, g_d, t_d) * 1e15, "--", color="gray",
            label=f"coax, $\\varepsilon_{{eff}}=\\varepsilon_0(1+\\varepsilon_r)/2$")
    ax.set_xlabel(xlabel); ax.set_ylabel("C (fF)")
    ax.set_title(f"{sw} sweep"); ax.grid(alpha=0.3); ax.legend(fontsize=8)

fig.suptitle("Single circular island: island$\\to$ground capacitance vs coaxial-edge model")
fig.savefig(f"{OUT}/island_cap_fits.png", dpi=130)
print(f"wrote {OUT}/island_cap_fits.png")

with open(f"{OUT}/island_cap_fit_report.json", "w") as f:
    json.dump({"eps_r": eps_r, "eps_eff_coplanar_over_eps0": eps_eff_coplanar / EPS0,
               "per_sweep": fit_report}, f, indent=2)
print(json.dumps(fit_report, indent=2))
