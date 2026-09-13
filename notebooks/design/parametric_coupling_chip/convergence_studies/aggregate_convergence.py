#!/usr/bin/env python
"""Convergence/reliability check for run_qc_sim.py's fitted (fr, Qc) at the
cs14_cl300 reference point: does the answer actually agree as solver_tol,
solver_maxits, mesh resolution, and AdaptiveMaxSamples are bracketed both
looser and tighter than our current defaults? See submit_convergence.sh's
header for the full rationale.

Run ON THE CLUSTER (per this repo's convention: generate sweep plots on the
cluster, not by pulling raw data back for local plotting) right after
submit_convergence.sh's array finishes:

    python aggregate_convergence.py
"""
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo/data root, one level up from this subfolder


def main():
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE, "out_conv", "results", "*.json"))):
        with open(f) as fh:
            rows.append(json.load(fh))
    if not rows:
        raise SystemExit("no results in out_conv/results -- has the sweep finished?")

    df = pd.DataFrame(rows)
    order = ["baseline", "tol_tighter", "tol_1e-5", "tol_1e-4", "tol_1e-4_maxits50",
            "tol_1e-3", "mesh_finer", "mesh_1.5x", "mesh_2x", "mesh_3x",
            "ams_more", "ams_200", "ams_120", "ams_80"]
    df["tag"] = pd.Categorical(df["tag"], categories=[t for t in order if t in df["tag"].values],
                               ordered=True)
    df = df.sort_values("tag")

    print(df[["tag", "solver_tol", "solver_maxits", "mesh_factor", "adaptive_max_samples",
              "fr_ghz", "Qr", "Qc", "Qi", "dip_depth_db", "solve_s"]].to_string(index=False))

    base = df[df["tag"] == "baseline"].iloc[0]
    print(f"\nbaseline: fr={base['fr_ghz']} GHz, Qc={base['Qc']:.1f}, solve_s={base['solve_s']:.1f}")

    df["fr_dev_pct"] = 100 * (df["fr_ghz"] - base["fr_ghz"]) / base["fr_ghz"]
    df["Qc_dev_pct"] = 100 * (df["Qc"] - base["Qc"]) / base["Qc"]
    print("\ndeviation from baseline:")
    print(df[["tag", "fr_dev_pct", "Qc_dev_pct", "solve_s"]].to_string(index=False))

    flagged = df[(df["fr_dev_pct"].abs() > 1.0) | (df["Qc_dev_pct"].abs() > 10.0)]
    if len(flagged):
        print(f"\nWARNING: {len(flagged)} point(s) deviate from baseline by >1% (fr) or "
             f">10% (Qc) -- these settings are NOT reliable / not converged:")
        print(flagged[["tag", "fr_dev_pct", "Qc_dev_pct"]].to_string(index=False))
    else:
        print("\nAll points agree with baseline within 1% (fr) / 10% (Qc) -- "
             "the fit is stable across this whole settings bracket.")

    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import thesis_style as ts
    ts.use()

    fig, axs = plt.subplots(2, 1, figsize=ts.figsize(1, ratio=1.2))
    x = np.arange(len(df))
    labels = df["tag"].astype(str).tolist()

    axs[0].axhline(base["fr_ghz"], ls="--", lw=0.8, color="gray")
    axs[0].plot(x, df["fr_ghz"], "o", color=ts.CYCLE[0])
    axs[0].set_ylabel(r"$f_r$ (GHz)")

    axs[1].axhline(base["Qc"], ls="--", lw=0.8, color="gray")
    axs[1].semilogy(x, df["Qc"], "o", color=ts.CYCLE[1])
    axs[1].set_ylabel(r"$Q_c$")
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(labels, rotation=60, ha="right", fontsize=6)

    ts.finalize(axs)
    ts.panel_labels(axs)
    IMAGES_DIR = os.path.join(HERE, "images")
    os.makedirs(IMAGES_DIR, exist_ok=True)
    ts.save(fig, "convergence_cs14_cl300", directory=IMAGES_DIR)
    print(f"\nwrote convergence_cs14_cl300.{{png,pdf}} in {HERE}")


if __name__ == "__main__":
    main()
