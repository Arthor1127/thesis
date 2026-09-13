#!/usr/bin/env python
"""Aggregates out_qc_eig/results/*.json (the eigenmode Qc sweep) into one CSV
and a thesis_style Qc(coupling_length) plot per coupling_space, plus a
side-by-side comparison against the driven-S21-fit reference (qc_prod_refit.csv,
where available) to show the eigenmode/driven offset across the whole grid,
not just the 3 points checked by hand this session.

Run ON THE CLUSTER: qsub submit_aggregate_eigenmode.sh
"""
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo/data root, one level up from this subfolder


def main():
    rows = []
    for jf in sorted(glob.glob(os.path.join(HERE, "out_qc_eig", "results", "*.json"))):
        with open(jf) as f:
            r = json.load(f)
        if not r.get("tag", "").endswith(("_eig", "_1port", "_t8")):
            continue
        rows.append(r)
    df = pd.DataFrame(rows)
    # Exclude cs14_cl300(_eig/_1port): a one-off validation/reference point tied
    # to the ORIGINAL driven-sim test geometry, not a member of the actual
    # production grid (coupling_space_um in {10,15,20,22.5,25,30}). Leaving it
    # in previously added a 7th category to a 6-color categorical palette,
    # silently making two real grid values (cs=10 and cs=30) share a color.
    GRID_SPACES = {10.0, 15.0, 20.0, 22.5, 25.0, 30.0}
    df = df[df["coupling_space_um"].isin(GRID_SPACES)].copy()
    df = df[["tag", "coupling_space_um", "coupling_length_um", "fr_ghz", "loaded_Q",
             "kappa_hz", "Qc"]].copy()
    df = df.sort_values(["coupling_space_um", "coupling_length_um"])
    out_csv = os.path.join(HERE, "qc_eigenmode.csv")
    df.to_csv(out_csv, index=False)
    print(f"wrote {out_csv} ({len(df)} points)")
    print(df.to_string(index=False))

    # Cross-check against the driven-S21-fit reference, where available
    # (qc_prod_refit.csv, from refit_qc_prod.py earlier this session).
    drv_path = os.path.join(HERE, "qc_prod_refit.csv")
    merged = None
    if os.path.exists(drv_path):
        drv = pd.read_csv(drv_path)
        drv = drv[(~drv["edge_flagged"].fillna(False).astype(bool)) &
                  (~drv["weak_local_dip"].fillna(False).astype(bool)) &
                  (drv["source"] == "out_qc_prod")].copy()
        drv = drv[["coupling_space_um", "coupling_length_um", "fr_ghz_refit", "Qc_refit"]]
        merged = df.merge(drv, on=["coupling_space_um", "coupling_length_um"], how="inner",
                          suffixes=("_eig", "_drv"))
        if len(merged):
            merged["Qc_ratio_eig_over_drv"] = merged["Qc"] / merged["Qc_refit"]
            merged["fr_diff_mhz"] = (merged["fr_ghz"] - merged["fr_ghz_refit"]) * 1e3
            out_cmp = os.path.join(HERE, "qc_eigenmode_vs_driven.csv")
            merged.to_csv(out_cmp, index=False)
            print(f"\nwrote {out_cmp} ({len(merged)} points with both eigenmode AND "
                  f"trusted driven-fit Qc)")
            print(merged[["tag_eig" if "tag_eig" in merged.columns else "tag",
                          "coupling_space_um", "coupling_length_um", "Qc", "Qc_refit",
                          "Qc_ratio_eig_over_drv", "fr_diff_mhz"]].to_string(index=False))
            ratios = merged["Qc_ratio_eig_over_drv"]
            print(f"\nQc_ratio (eigenmode/driven): mean={ratios.mean():.3f} "
                  f"std={ratios.std():.3f} min={ratios.min():.3f} max={ratios.max():.3f} "
                  f"n={len(ratios)}")
    else:
        print(f"\n{drv_path} not found -- skipping driven-fit cross-check")

    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import thesis_style as ts
    ts.use()

    fig, ax = plt.subplots(figsize=ts.figsize(1))
    for cs, grp in df.groupby("coupling_space_um"):
        grp = grp.sort_values("coupling_length_um")
        ax.semilogy(grp["coupling_length_um"], grp["Qc"], "o-", label=fr"$s={cs:g}\,\mu$m")
    ax.set_xlabel(r"coupling length $\ell$ ($\mu$m)")
    ax.set_ylabel(r"$Q_c$ (eigenmode)")
    ax.legend(fontsize=6, ncol=2)
    ts.finalize([ax])
    IMAGES_DIR = os.path.join(HERE, "images")
    os.makedirs(IMAGES_DIR, exist_ok=True)
    ts.save(fig, "qc_eigenmode", directory=IMAGES_DIR)
    print(f"\nwrote qc_eigenmode.{{png,pdf}}")

    if merged is not None and len(merged):
        fig2, axs = plt.subplots(2, 1, figsize=ts.figsize(1, ratio=1.2), sharex=True)
        for cs, grp in merged.groupby("coupling_space_um"):
            grp = grp.sort_values("coupling_length_um")
            axs[0].semilogy(grp["coupling_length_um"], grp["Qc"], "o-",
                            label=fr"$s={cs:g}\,\mu$m (eig)")
            axs[0].semilogy(grp["coupling_length_um"], grp["Qc_refit"], "s--",
                            color=axs[0].lines[-1].get_color(), alpha=0.5)
        axs[0].set_ylabel(r"$Q_c$")
        axs[0].legend(fontsize=6, ncol=2)
        for cs, grp in merged.groupby("coupling_space_um"):
            grp = grp.sort_values("coupling_length_um")
            axs[1].plot(grp["coupling_length_um"], grp["Qc_ratio_eig_over_drv"], "o-")
        axs[1].axhline(1.0, ls="--", lw=0.8, color="gray")
        axs[1].set_xlabel(r"coupling length $\ell$ ($\mu$m)")
        axs[1].set_ylabel(r"$Q_c^{eig}/Q_c^{drv}$")
        ts.finalize(axs)
        ts.panel_labels(axs)
        ts.save(fig2, "qc_eigenmode_vs_driven", directory=IMAGES_DIR)
        print(f"wrote qc_eigenmode_vs_driven.{{png,pdf}}")


if __name__ == "__main__":
    main()
