#!/usr/bin/env python
"""Validate DesignHelper's finite-metal-thickness correction (Eqs. B6-B7 of
compute_cpw_geometric_parameters) against FULL3D Palace sims at a few
thicknesses, for each representative (width, gap) point in the thickness
sweep (gen_thickness_params.py, submit_*_thickness.sh).

For each (width, gap, thickness) triple, pairs the capacitance and
inductance results, computes Z0_sim/v_sim exactly as aggregate_and_plot.py
does, and compares against ``build_design.analytic_reference(width, gap,
thickness_nm=thickness)`` evaluated AT THAT SAME THICKNESS -- unlike the main
sweep (which fixes REF_THICKNESS_NM=100nm on the analytic side while Palace
uses zero-thickness metal), this is an apples-to-apples thickness-resolved
comparison in both directions.

Runs either locally or on the cluster (thesis_style is installed in both, see
aggregate_and_plot.py's docstring).
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

from build_design import analytic_reference

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(results_dir):
    rows = []
    for f in sorted(glob.glob(os.path.join(results_dir, "*.json"))):
        with open(f) as fh:
            rows.append(json.load(fh))
    return rows


def build_table(cap_dir, ind_dir):
    cap_rows = {r["tag"]: r for r in _load(cap_dir)}
    ind_rows = {r["tag"]: r for r in _load(ind_dir)}
    tags = sorted(set(cap_rows) & set(ind_rows))
    missing = sorted(set(cap_rows) ^ set(ind_rows))
    if missing:
        print(f"WARNING: {len(missing)} tag(s) present in only one sweep, skipped: {missing}")

    recs = []
    for tag in tags:
        c, l = cap_rows[tag], ind_rows[tag]
        assert c["full3d"] and l["full3d"], f"{tag}: expected FULL3D results from the thickness sweep"
        w, g = c["cpw_width_um"], c["cpw_gap_um"]
        t_nm = c["metal_thickness_nm"]
        assert abs(t_nm - l["metal_thickness_nm"]) < 1e-6, f"{tag}: thickness mismatch cap vs ind"

        Cprime, Lprime = c["Cprime_F_per_m"], l["Lprime_H_per_m"]
        Z0_sim = np.sqrt(Lprime / Cprime)
        v_sim = 1.0 / np.sqrt(Lprime * Cprime)

        ref = analytic_reference(w, g, thickness_nm=t_nm)
        Z0_ref = ref["Z0"].to("ohm").magnitude
        v_ref = ref["v"].to("m/s").magnitude

        recs.append(dict(
            tag=tag, cpw_width_um=w, cpw_gap_um=g, thickness_nm=t_nm,
            Z0_sim_ohm=Z0_sim, v_sim_m_per_s=v_sim, Z0_ref_ohm=Z0_ref, v_ref_m_per_s=v_ref,
            Z0_rel_err_pct=100.0 * (Z0_sim - Z0_ref) / Z0_ref,
            v_rel_err_pct=100.0 * (v_sim - v_ref) / v_ref,
        ))
    return pd.DataFrame.from_records(recs).sort_values(["cpw_width_um", "cpw_gap_um", "thickness_nm"])


def make_plots(df, outdir):
    import matplotlib.pyplot as plt
    import thesis_style as ts
    ts.use()

    fig, axs = plt.subplots(1, 2, figsize=ts.figsize(2, ratio=0.5))
    points = sorted(df[["cpw_width_um", "cpw_gap_um"]].drop_duplicates().itertuples(index=False))
    for ax, (col_sim, col_ref, ylabel) in zip(
        axs, [("Z0_sim_ohm", "Z0_ref_ohm", r"$Z_0$ (Ohm)"),
              ("v_sim_m_per_s", "v_ref_m_per_s", r"$v$ (m/s)")]):
        for i, (w, g) in enumerate(points):
            sub = df[(df["cpw_width_um"] == w) & (df["cpw_gap_um"] == g)].sort_values("thickness_nm")
            color = ts.CYCLE[i % len(ts.CYCLE)]
            ax.plot(sub["thickness_nm"], sub[col_ref], "-", color=color,
                   label=f"w={w:g} g={g:g} um (model)")
            ax.plot(sub["thickness_nm"], sub[col_sim], "o", color=color,
                   label=f"w={w:g} g={g:g} um (Palace)")
        ax.set_xlabel(r"film thickness $t$ (nm)")
        ax.set_ylabel(ylabel)
    axs[0].legend(fontsize=6)
    ts.finalize(axs)
    ts.panel_labels(axs)
    ts.save(fig, "cpw_thickness_correction", directory=outdir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap-dir", default=os.path.join(HERE, "out_cap_thickness", "results"))
    ap.add_argument("--ind-dir", default=os.path.join(HERE, "out_ind_thickness", "results"))
    ap.add_argument("--out-csv", default=os.path.join(HERE, "cpw_thickness_results.csv"))
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    df = build_table(args.cap_dir, args.ind_dir)
    df.to_csv(args.out_csv, index=False)
    print(f"wrote {args.out_csv} ({len(df)} points)")
    print(df.to_string(index=False))

    if not args.no_plots:
        make_plots(df, HERE)


if __name__ == "__main__":
    main()
