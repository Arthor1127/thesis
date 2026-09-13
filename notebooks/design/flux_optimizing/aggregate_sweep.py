#!/usr/bin/env python
"""Combine sweep/shard_*.csv -> sweep_results.csv, print the best configs, and
render pairwise 2-D colormaps (max-projected over the other three axes).

    python aggregate_sweep.py --sweep-dir sweep
"""
import argparse
import glob
import itertools
import os

import numpy as np
import pandas as pd

AXES = ['route_clearance_um', 'segment_length_um', 'squid_gap_um',
        'arm_semi_length_um', 'arm_width_um']
NICE = {'route_clearance_um': 'route clearance (um)',
        'segment_length_um': 'inter-loop lead (um)',
        'squid_gap_um': 'loop height (um)',
        'arm_semi_length_um': 'loop half-width (um)',
        'arm_width_um': 'lead width (um)'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sweep-dir', default='sweep')
    ap.add_argument('--metric', default='Phi0_per_mA')
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.sweep_dir, 'shard_*.csv')))
    if not files:
        raise SystemExit(f"no shard_*.csv in {args.sweep_dir}")
    parts = []
    for f in files:
        try:
            parts.append(pd.read_csv(f, on_bad_lines='skip'))
        except Exception as e:
            print(f"  (skipped unreadable {f}: {e})")
    df = pd.concat(parts, ignore_index=True)
    df = df.dropna(subset=AXES)
    df = df.drop_duplicates(subset=AXES, keep='last').sort_values(AXES).reset_index(drop=True)

    import sweep_flux
    n_total = len(sweep_flux.all_points())
    ok = df[args.metric].notna() & (df['error'].isna() | (df['error'].astype(str) == ''))
    print(f"{len(df)}/{n_total} points present, {ok.sum()} succeeded, "
          f"{(~ok).sum()} failed/NaN")
    if 'converged' in df:
        nconv = (df.loc[ok, 'converged'].astype(str).str.lower() == 'false').sum()
        if nconv:
            print(f"  ({nconv} succeeded but solver flagged not-converged -- inspect)")

    out_csv = os.path.join(args.sweep_dir, 'sweep_results.csv')
    df.to_csv(out_csv, index=False)
    print(f"wrote {out_csv}")

    g = df[ok].copy()
    print(f"\ntop 15 by {args.metric}:")
    cols = AXES + [args.metric, 'flux_per_amp', 'n_ksp', 'seconds']
    print(g.nlargest(15, args.metric)[cols].to_string(index=False))
    print(f"\nbottom 5:")
    print(g.nsmallest(5, args.metric)[cols].to_string(index=False))

    # 1-D marginal maxima (best achievable at each value of one axis)
    print("\n1-D marginal max of", args.metric, "(best over all other axes):")
    for a in AXES:
        s = g.groupby(a)[args.metric].max()
        print(f"  {NICE[a]:24s} " + "  ".join(f"{v:.3g}@{k:g}" for k, v in s.items()))

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import thesis_style as ts
        import thesis_cmaps  # noqa: F401  -- registers 'thesis_cmap'
        ts.use()
    except Exception:
        print("\n(matplotlib/thesis_style unavailable -- skipping colormaps)")
        return

    pairs = list(itertools.combinations(AXES, 2))
    ncol = 5
    nrow = int(np.ceil(len(pairs) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.6 * nrow), squeeze=False)
    vmax = g[args.metric].max()
    vmin = g[args.metric].min()
    for ax, (ax1, ax2) in zip(axs.ravel(), pairs):
        piv = g.pivot_table(index=ax2, columns=ax1, values=args.metric, aggfunc='max')
        im = ax.imshow(piv.values, origin='lower', aspect='auto', cmap='thesis_cmap',
                       vmin=vmin, vmax=vmax,
                       extent=[piv.columns.min(), piv.columns.max(),
                               piv.index.min(), piv.index.max()])
        ax.set_xlabel(NICE[ax1]); ax.set_ylabel(NICE[ax2])
        fig.colorbar(im, ax=ax, shrink=0.85)
    for ax in axs.ravel()[len(pairs):]:
        ax.axis('off')
    ts.finalize(axs.ravel()[:len(pairs)])
    fig.suptitle(f"max {args.metric} projected onto each axis pair", y=1.0)
    fig.tight_layout()
    png = os.path.join(args.sweep_dir, 'sweep_colormaps.png')
    fig.savefig(png, dpi=110, bbox_inches='tight', facecolor='white')
    print(f"\nwrote {png}")


if __name__ == '__main__':
    main()
