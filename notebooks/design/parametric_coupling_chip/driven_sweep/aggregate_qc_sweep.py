#!/usr/bin/env python
"""Combine the Qc sweep results (run_qc_sim.py via submit_qc_sweep.sh) into
f_res(coupling_space_um, coupling_length_um) and Qc(coupling_space_um,
coupling_length_um) plots. Uses thesis_style/thesis_cmaps
(packages/qcg-quantum-design, per directives/figure-style.md), which are
installed in this machine's own conda envs AND (as of this session) the
cluster's qcg-quantum-design env -- so this can run either locally after
pulling out_qc/results/ back, or directly on the cluster right after the
sweep finishes.

    python aggregate_qc_sweep.py
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo/data root, one level up from this subfolder


def _load(results_dir):
    rows = []
    for f in sorted(glob.glob(os.path.join(results_dir, "*.json"))):
        with open(f) as fh:
            rows.append(json.load(fh))
    return rows


def make_plots(df, outdir):
    import matplotlib.pyplot as plt
    import thesis_style as ts
    ts.use()

    spaces = sorted(df["coupling_space_um"].unique())
    fig, axs = plt.subplots(1, 2, figsize=ts.figsize(2, ratio=0.5))
    for i, cs in enumerate(spaces):
        sub = df[df["coupling_space_um"] == cs].sort_values("coupling_length_um")
        color = ts.CYCLE[i % len(ts.CYCLE)]
        axs[0].plot(sub["coupling_length_um"], sub["fr_ghz"], "o-", color=color,
                   label=f"space={cs:g} um")
        axs[1].semilogy(sub["coupling_length_um"], sub["Qc"], "o-", color=color,
                       label=f"space={cs:g} um")
    axs[0].set_xlabel(r"coupling length (um)")
    axs[0].set_ylabel(r"$f_r$ (GHz)")
    axs[1].set_xlabel(r"coupling length (um)")
    axs[1].set_ylabel(r"$Q_c$")
    axs[0].legend(fontsize=6)
    ts.finalize(axs)
    ts.panel_labels(axs)
    ts.save(fig, "qc_sweep_vs_length", directory=outdir)

    # Qc heatmap over the full (space, length) grid
    lengths = sorted(df["coupling_length_um"].unique())
    grid = df.pivot(index="coupling_space_um", columns="coupling_length_um", values="Qc") \
            .reindex(index=spaces, columns=lengths)
    fig2, ax2 = plt.subplots(figsize=ts.figsize(1, ratio=0.8))
    im = ax2.imshow(np.log10(grid.values), origin="lower", aspect="auto", cmap="thesis_cmap",
                    extent=[min(lengths), max(lengths), min(spaces), max(spaces)])
    ax2.set_xlabel(r"coupling length (um)")
    ax2.set_ylabel(r"coupling space (um)")
    cbar = fig2.colorbar(im, ax=ax2)
    cbar.set_label(r"$\log_{10} Q_c$")
    ts.save(fig2, "qc_sweep_heatmap", directory=outdir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default=os.path.join(HERE, "out_qc", "results"))
    ap.add_argument("--out-csv", default=os.path.join(HERE, "qc_sweep_results.csv"))
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    rows = _load(args.results_dir)
    if not rows:
        raise SystemExit(f"no results in {args.results_dir} -- pull them back from the "
                         "cluster first, or run this on the cluster directly")
    df = pd.DataFrame(rows).sort_values(["coupling_space_um", "coupling_length_um"])
    df.to_csv(args.out_csv, index=False)
    print(f"wrote {args.out_csv} ({len(df)} points)")
    print(df[["coupling_space_um", "coupling_length_um", "fr_ghz", "Qr", "Qc", "Qi",
             "dip_freq_ghz", "dip_depth_db"]].to_string(index=False))

    shallow = df[df["dip_depth_db"] > -3.0]
    if len(shallow):
        print(f"\nWARNING: {len(shallow)} point(s) have dip_depth_db > -3dB -- the "
             "frequency window likely missed the true resonance there (see "
             "run_qc_sim.py's FREQUENCY WINDOW note); Qr/Qc/Qi for these are not "
             "trustworthy. Check the corresponding *_fit.png:")
        print(shallow[["coupling_space_um", "coupling_length_um", "dip_depth_db"]]
             .to_string(index=False))

    if not args.no_plots:
        make_plots(df, HERE)


if __name__ == "__main__":
    main()
