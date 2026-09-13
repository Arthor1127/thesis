#!/usr/bin/env python
"""Collect submit_feedline_sweep.sh's eigenmode results and separate the
resonator's own mode from the feedline's cavity harmonics.

The discriminator needs no field analysis: an open-open feedline cavity has
harmonics at f_n = n*v/(2L), so every feedline mode moves as 1/L across the
sweep, while the resonator's own mode (fixed length, fixed coupling) sits at
essentially the same frequency at every point. Plot all modes vs feedline
length and the ~horizontal line is the resonator; the fanning family is the
feedline.

Uses thesis_style/thesis_cmaps per directives/figure-style.md.

    python aggregate_feedline_sweep.py                 # after pulling out_fdln/ back
    python aggregate_feedline_sweep.py --no-plots      # table only
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo/data root, one level up from this subfolder


def load(results_dir):
    rows = []
    for f in sorted(glob.glob(os.path.join(results_dir, "*", "resonance.json"))):
        with open(f) as fh:
            r = json.load(fh)
        r["_path"] = f
        rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default=os.path.join(HERE, "out_fdln"))
    ap.add_argument("--out-csv", default=os.path.join(HERE, "feedline_sweep_modes.csv"))
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--resonator-guess-ghz", type=float, default=8.29,
                    help="Anchor for picking the resonator mode at each length. The "
                        "12mm point (job 1938824) had its resonator-localised mode -- "
                        "confirmed from the E-field plots -- at 8.291GHz.")
    args = ap.parse_args()

    rows = load(args.results_dir)
    if not rows:
        raise SystemExit(f"no resonance.json under {args.results_dir}")

    recs = []
    for r in rows:
        if r.get("status") != "success":
            print(f"WARNING: {r['_path']} status={r.get('status')} -- skipped")
            continue
        L = r["feedline_length_mm"]
        for i, (f, q) in enumerate(zip(r["all_modes_ghz"], r["all_q_factors"])):
            recs.append({"feedline_mm": L, "half_span_um": r["chip_half_span_um"],
                         "mode_index": i, "f_ghz": f, "Q": q})
    df = pd.DataFrame(recs).sort_values(["feedline_mm", "f_ghz"])
    df.to_csv(args.out_csv, index=False)
    print(f"wrote {args.out_csv} ({len(df)} modes across "
          f"{df['feedline_mm'].nunique()} feedline lengths)\n")

    # per-length: the mode closest to the resonator anchor, and its nearest neighbour
    print(f"{'feedline':>9} {'f_res*':>8} {'nearest other':>14} {'gap':>8}   all modes (GHz)")
    summary = []
    for L, sub in df.groupby("feedline_mm"):
        fs = np.sort(sub["f_ghz"].values)
        i = int(np.argmin(np.abs(fs - args.resonator_guess_ghz)))
        f_res = fs[i]
        others = np.delete(fs, i)
        if len(others):
            j = int(np.argmin(np.abs(others - f_res)))
            gap = abs(others[j] - f_res)
            near = others[j]
        else:
            gap, near = np.inf, np.nan
        summary.append({"feedline_mm": L, "f_res_ghz": f_res,
                        "nearest_other_ghz": near, "gap_ghz": gap})
        print(f"{L:8.1f}mm {f_res:8.3f} {near:14.3f} {gap:8.3f}   "
              + ", ".join(f"{f:.3f}" for f in fs))
    sdf = pd.DataFrame(summary)
    best = sdf.loc[sdf["gap_ghz"].idxmax()]
    print(f"\n* 'f_res' = mode nearest {args.resonator_guess_ghz}GHz; verify against the "
          f"eig*_ErealMag.png field plots before trusting it.")
    print(f"Largest resonator-to-neighbour gap: {best['gap_ghz']:.3f}GHz at "
          f"feedline = {best['feedline_mm']:.1f}mm")

    if args.no_plots:
        return

    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import thesis_style as ts
    ts.use()

    fig, axs = plt.subplots(1, 2, figsize=ts.figsize(2, ratio=0.5))
    axs[0].plot(df["feedline_mm"], df["f_ghz"], "o", color=ts.CYCLE[0], label="eigenmodes")
    axs[0].axhline(args.resonator_guess_ghz, ls="--", lw=0.8, color=ts.CYCLE[1],
                   label=r"resonator anchor")
    # overlay the expected 1/L feedline harmonic family for reference
    Lmm = np.linspace(df["feedline_mm"].min() * 0.9, df["feedline_mm"].max() * 1.1, 200)
    v = 1.207e8  # m/s, same phase velocity the resonator length targets 8GHz with
    for n in range(1, 8):
        f_n = n * v / (2 * Lmm * 1e-3) / 1e9
        if f_n.max() < df["f_ghz"].min() or f_n.min() > df["f_ghz"].max():
            continue
        axs[0].plot(Lmm, f_n, "-", lw=0.7, color=ts.CYCLE[2], alpha=0.6,
                    label=r"$nv/2L$" if n == 1 else None)
    axs[0].set_ylim(df["f_ghz"].min() - 0.5, df["f_ghz"].max() + 0.5)
    axs[0].set_xlabel("feedline length (mm)")
    axs[0].set_ylabel(r"$f$ (GHz)")
    axs[0].legend(fontsize=6)

    axs[1].plot(sdf["feedline_mm"], sdf["gap_ghz"], "o-", color=ts.CYCLE[0])
    axs[1].set_xlabel("feedline length (mm)")
    axs[1].set_ylabel("resonator-to-nearest-mode gap (GHz)")

    ts.finalize(axs)
    ts.panel_labels(axs)
    IMAGES_DIR = os.path.join(HERE, "images")
    os.makedirs(IMAGES_DIR, exist_ok=True)
    ts.save(fig, "feedline_sweep_modes", directory=IMAGES_DIR)
    print(f"\nwrote feedline_sweep_modes.{{png,pdf}} in {HERE}")


if __name__ == "__main__":
    main()
