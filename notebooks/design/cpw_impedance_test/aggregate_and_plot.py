#!/usr/bin/env python
"""Combine the capacitance + self-inductance sweep results into Z0(cpw_width,
cpw_gap) and v(cpw_width, cpw_gap), and compare against
``parametric_coupling.design_helper.DesignHelper.compute_cpw_geometric_parameters``.

Needs ``thesis_style``/``thesis_cmaps`` (packages/qcg-quantum-design), which
are pip-installed into this machine's own `quantum` / `qcg-quantum-design`
conda envs AND (as of `submit_install_thesis_style.sh`) the cluster's
`qcg-quantum-design` env -- so this can run either locally after pulling
results back, e.g.:

    rsync -avz arturo.ruiz@<cluster>:circuit_design/cpw_impedance_test/out_cap/results/ out_cap/results/
    rsync -avz arturo.ruiz@<cluster>:circuit_design/cpw_impedance_test/out_ind/results/ out_ind/results/
    python aggregate_and_plot.py   # in the qcg-quantum-design (or quantum) conda env

or directly on the cluster, in ~/circuit_design/cpw_impedance_test with the
qcg-quantum-design env activated, right after the sweep finishes.
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
        w, g = c["cpw_width_um"], c["cpw_gap_um"]
        Cprime = c["Cprime_F_per_m"]           # F/m
        Lprime = l["Lprime_H_per_m"]           # H/m
        Z0_sim = np.sqrt(Lprime / Cprime)
        v_sim = 1.0 / np.sqrt(Lprime * Cprime)

        ref = analytic_reference(w, g)
        Z0_ref = ref["Z0"].to("ohm").magnitude
        v_ref = ref["v"].to("m/s").magnitude

        recs.append(dict(
            tag=tag, cpw_width_um=w, cpw_gap_um=g,
            Cprime_fF_per_um=c["Cprime_fF_per_um"], Lprime_nH_per_um=l["Lprime_nH_per_um"],
            Z0_sim_ohm=Z0_sim, v_sim_m_per_s=v_sim,
            Z0_ref_ohm=Z0_ref, v_ref_m_per_s=v_ref,
            Z0_rel_err_pct=100.0 * (Z0_sim - Z0_ref) / Z0_ref,
            v_rel_err_pct=100.0 * (v_sim - v_ref) / v_ref,
        ))
    return pd.DataFrame.from_records(recs).sort_values(["cpw_gap_um", "cpw_width_um"])


def make_plots(df, outdir):
    import matplotlib.pyplot as plt
    import thesis_style as ts
    ts.use()

    fig, axs = plt.subplots(1, 2, figsize=ts.figsize(2, ratio=0.5))

    gaps = sorted(df["cpw_gap_um"].unique())
    for ax, (col_sim, col_ref, ylabel) in zip(
        axs, [("Z0_sim_ohm", "Z0_ref_ohm", r"$Z_0$ (Ohm)"),
              ("v_sim_m_per_s", "v_ref_m_per_s", r"$v$ (m/s)")]):
        for i, g in enumerate(gaps):
            sub = df[df["cpw_gap_um"] == g].sort_values("cpw_width_um")
            color = ts.CYCLE[i % len(ts.CYCLE)]
            ax.plot(sub["cpw_width_um"], sub[col_ref], "-", color=color,
                   label=f"g={g:g} um (model)")
            ax.plot(sub["cpw_width_um"], sub[col_sim], "o", color=color,
                   label=f"g={g:g} um (Palace)")
        ax.set_xlabel(r"$W$ (um)")
        ax.set_ylabel(ylabel)

    axs[0].legend(fontsize=6, ncol=2)
    ts.finalize(axs)
    ts.panel_labels(axs)
    ts.save(fig, "cpw_impedance_validation", directory=outdir)

    # relative-error heatmaps
    fig2, axs2 = plt.subplots(1, 2, figsize=ts.figsize(2, ratio=0.5))
    widths = sorted(df["cpw_width_um"].unique())
    for ax, col, title in zip(axs2, ["Z0_rel_err_pct", "v_rel_err_pct"],
                              [r"$Z_0$ rel. error (%)", r"$v$ rel. error (%)"]):
        grid = df.pivot(index="cpw_gap_um", columns="cpw_width_um", values=col) \
                .reindex(index=gaps, columns=widths)
        im = ax.imshow(grid.values, origin="lower", aspect="auto", cmap="thesis_cmap_diverging",
                       vmin=-np.nanmax(np.abs(grid.values)), vmax=np.nanmax(np.abs(grid.values)),
                       extent=[min(widths), max(widths), min(gaps), max(gaps)])
        ax.set_xlabel(r"$W$ (um)")
        ax.set_ylabel(r"$G$ (um)")
        ax.set_title(title)
        fig2.colorbar(im, ax=ax)
    ts.save(fig2, "cpw_impedance_rel_error", directory=outdir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap-dir", default=os.path.join(HERE, "out_cap", "results"))
    ap.add_argument("--ind-dir", default=os.path.join(HERE, "out_ind", "results"))
    ap.add_argument("--out-csv", default=os.path.join(HERE, "cpw_impedance_results.csv"))
    ap.add_argument("--no-plots", action="store_true",
                    help="only write the CSV (skip thesis_style, e.g. to sanity-check "
                         "on the cluster before pulling data back)")
    args = ap.parse_args()

    df = build_table(args.cap_dir, args.ind_dir)
    df.to_csv(args.out_csv, index=False)
    print(f"wrote {args.out_csv} ({len(df)} points)")
    print(df[["cpw_width_um", "cpw_gap_um", "Z0_sim_ohm", "Z0_ref_ohm", "Z0_rel_err_pct",
             "v_rel_err_pct"]].to_string(index=False))

    if not args.no_plots:
        make_plots(df, HERE)


if __name__ == "__main__":
    main()
