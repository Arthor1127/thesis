#!/usr/bin/env python
"""Tests the mode-hybridization hypothesis quantitatively across the whole
production grid: for each point with both an eigenmode result (which saves
the full nearby-mode list) and a trusted driven-fit Qc, compute the gap to
the nearest competing eigenmode and correlate it against
Qc(eigenmode)/Qc(driven-fit). If mode-crowding is really what's corrupting
points like cs25_cl200 (mode gap 147MHz, the worst outlier at ratio 0.58),
smaller gaps should predict ratios further from 1.

Run ON THE CLUSTER: qsub submit_mode_gap_correlation.sh
"""
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo/data root, one level up from this subfolder


def main():
    gap_rows = []
    for jf in sorted(glob.glob(os.path.join(HERE, "out_qc_eig", "results", "*_eig.json"))):
        with open(jf) as f:
            r = json.load(f)
        modes = r.get("modes", [])
        if len(modes) < 2:
            continue
        fr = r["fr_ghz"]
        # nearest OTHER mode (by frequency), not necessarily modes[1]
        others = [m["f_ghz"] for m in modes if abs(m["f_ghz"] - fr) > 1e-6]
        if not others:
            continue
        nearest_gap_mhz = min(abs(f - fr) for f in others) * 1e3
        gap_rows.append({"tag": r["tag"], "coupling_space_um": r["coupling_space_um"],
                         "coupling_length_um": r["coupling_length_um"],
                         "fr_ghz": fr, "nearest_mode_gap_mhz": nearest_gap_mhz})
    gaps = pd.DataFrame(gap_rows)
    print(f"mode-gap data for {len(gaps)} points")

    cmp_path = os.path.join(HERE, "qc_eigenmode_vs_driven.csv")
    if not os.path.exists(cmp_path):
        print(f"{cmp_path} not found -- writing gap table only")
        gaps.to_csv(os.path.join(HERE, "mode_gaps.csv"), index=False)
        return

    cmp = pd.read_csv(cmp_path)
    merged = cmp.merge(gaps, on=["coupling_space_um", "coupling_length_um"], how="inner",
                       suffixes=("", "_gapcalc"))
    merged = merged.sort_values("nearest_mode_gap_mhz")
    out_csv = os.path.join(HERE, "mode_gap_vs_ratio.csv")
    merged.to_csv(out_csv, index=False)
    print(f"wrote {out_csv} ({len(merged)} points)")
    print(merged[["tag", "nearest_mode_gap_mhz", "Qc_ratio_eig_over_drv"]].to_string(index=False))

    if len(merged) >= 3:
        corr = merged["nearest_mode_gap_mhz"].corr(merged["Qc_ratio_eig_over_drv"])
        abs_dev = (merged["Qc_ratio_eig_over_drv"] - 1).abs()
        corr_abs = merged["nearest_mode_gap_mhz"].corr(abs_dev)
        print(f"\nPearson corr(gap_mhz, ratio) = {corr:.3f}")
        print(f"Pearson corr(gap_mhz, |ratio-1|) = {corr_abs:.3f}  "
              f"(negative = smaller gap -> bigger deviation, as the hybridization "
              f"hypothesis predicts)")

    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import thesis_style as ts
    ts.use()

    fig, ax = plt.subplots(figsize=ts.figsize(1))
    ax.scatter(merged["nearest_mode_gap_mhz"], merged["Qc_ratio_eig_over_drv"],
              color=ts.CYCLE[0], zorder=3)
    for _, row in merged.iterrows():
        ax.annotate(row["tag"].replace("_eig", ""),
                   (row["nearest_mode_gap_mhz"], row["Qc_ratio_eig_over_drv"]),
                   fontsize=5, alpha=0.7, xytext=(3, 3), textcoords="offset points")
    ax.axhline(1.0, ls="--", lw=0.8, color="gray")
    ax.set_xlabel("nearest competing-mode gap (MHz)")
    ax.set_ylabel(r"$Q_c^{eig}/Q_c^{drv}$")
    ts.finalize([ax])
    IMAGES_DIR = os.path.join(HERE, "images")
    os.makedirs(IMAGES_DIR, exist_ok=True)
    ts.save(fig, "mode_gap_vs_ratio", directory=IMAGES_DIR)
    print("wrote mode_gap_vs_ratio.{png,pdf}")


if __name__ == "__main__":
    main()
