#!/usr/bin/env python
"""Verify that C_total(length) and L_total(length) are actually linear in
length -- the premise behind treating a single length's C_total/length_um (or
L_total/integration_length_um) as "the" per-unit-length C'/L'. Reads the
length-linearity sweep (gen_lengthcheck_params.py, submit_*_lengthcheck.sh).

For each (cpw_width_um, cpw_gap_um) point, fits a line through
(length, C_total) [capacitance] and (integration_length, L_total)
[inductance -- the flux integral is margin-trimmed, so that trimmed length,
not the physical trace length, is the fair linear variable] and reports:

  - the fit slope, i.e. C'/L' as measured from the length DEPENDENCE alone,
    independent of any single length's fringe;
  - how that compares to Cprime/Lprime as reported by each individual length
    (they should agree once fringe is genuinely a length-independent additive
    constant);
  - the intercept, as a fraction of a typical total -- large means the "single
    length, divide by length" shortcut used by the main sweep carries a
    non-negligible systematic at the lengths it uses.

Runs either locally or on the cluster (thesis_style is installed in both, see
aggregate_and_plot.py's docstring).
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(results_dir):
    rows = []
    for f in sorted(glob.glob(os.path.join(results_dir, "*.json"))):
        with open(f) as fh:
            rows.append(json.load(fh))
    return rows


def _fit_line(x, y):
    """slope, intercept, and max |residual|/max(|y|) for a straight-line fit."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    max_rel_resid = float(np.max(np.abs(resid)) / np.max(np.abs(y)))
    return float(slope), float(intercept), max_rel_resid


def analyze(cap_dir, ind_dir):
    cap_rows = _load(cap_dir)
    ind_rows = _load(ind_dir)
    if not cap_rows or not ind_rows:
        raise SystemExit(f"no results in {cap_dir} and/or {ind_dir} -- "
                         "pull them back from the cluster first (see README.md)")

    cap_df = pd.DataFrame(cap_rows)
    ind_df = pd.DataFrame(ind_rows)

    records = []
    for (w, g), grp in cap_df.groupby(["cpw_width_um", "cpw_gap_um"]):
        grp = grp.sort_values("length_um")
        slope_C, intercept_C, resid_C = _fit_line(grp["length_um"], grp["C_total_F"])
        # slope is F per um -> convert to fF/um for readability
        slope_C_fF_per_um = slope_C * 1e15
        typical_C = grp["C_total_F"].median()
        records.append(dict(
            quantity="C", cpw_width_um=w, cpw_gap_um=g, n_lengths=len(grp),
            slope_fF_per_um=slope_C_fF_per_um,
            intercept_fF=intercept_C * 1e15,
            intercept_frac_of_typical=abs(intercept_C) / typical_C,
            max_rel_fit_residual=resid_C,
            single_point_Cprime_fF_per_um=grp["Cprime_fF_per_um"].tolist(),
        ))

    for (w, g), grp in ind_df.groupby(["cpw_width_um", "cpw_gap_um"]):
        grp = grp.sort_values("integration_length_um")
        slope_L, intercept_L, resid_L = _fit_line(grp["integration_length_um"], grp["L_total_H"])
        slope_L_nH_per_um = slope_L * 1e9
        typical_L = grp["L_total_H"].median()
        records.append(dict(
            quantity="L", cpw_width_um=w, cpw_gap_um=g, n_lengths=len(grp),
            slope_fF_per_um=None, slope_nH_per_um=slope_L_nH_per_um,
            intercept_nH=intercept_L * 1e9,
            intercept_frac_of_typical=abs(intercept_L) / typical_L,
            max_rel_fit_residual=resid_L,
            single_point_Lprime_nH_per_um=grp["Lprime_nH_per_um"].tolist(),
        ))

    return pd.DataFrame.from_records(records), cap_df, ind_df


def make_plots(cap_df, ind_df, outdir):
    import matplotlib.pyplot as plt
    import thesis_style as ts
    ts.use()

    fig, axs = plt.subplots(1, 2, figsize=ts.figsize(2, ratio=0.5))
    for i, ((w, g), grp) in enumerate(cap_df.groupby(["cpw_width_um", "cpw_gap_um"])):
        grp = grp.sort_values("length_um")
        color = ts.CYCLE[i % len(ts.CYCLE)]
        axs[0].plot(grp["length_um"], grp["C_total_fF"], "o-", color=color,
                   label=f"w={w:g} g={g:g} um")
    axs[0].set_xlabel(r"length (um)")
    axs[0].set_ylabel(r"$C_\mathrm{total}$ (fF)")
    axs[0].legend(fontsize=6)

    for i, ((w, g), grp) in enumerate(ind_df.groupby(["cpw_width_um", "cpw_gap_um"])):
        grp = grp.sort_values("integration_length_um")
        color = ts.CYCLE[i % len(ts.CYCLE)]
        axs[1].plot(grp["integration_length_um"], grp["L_total_H"] * 1e9, "o-", color=color,
                   label=f"w={w:g} g={g:g} um")
    axs[1].set_xlabel(r"integration length (um)")
    axs[1].set_ylabel(r"$L_\mathrm{total}$ (nH)")
    axs[1].legend(fontsize=6)

    ts.finalize(axs)
    ts.panel_labels(axs)
    ts.save(fig, "cpw_length_linearity", directory=outdir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap-dir", default=os.path.join(HERE, "out_cap_lengthcheck", "results"))
    ap.add_argument("--ind-dir", default=os.path.join(HERE, "out_ind_lengthcheck", "results"))
    ap.add_argument("--out-csv", default=os.path.join(HERE, "cpw_length_linearity.csv"))
    ap.add_argument("--no-plots", action="store_true")
    args = ap.parse_args()

    summary, cap_df, ind_df = analyze(args.cap_dir, args.ind_dir)
    summary.to_csv(args.out_csv, index=False)
    print(f"wrote {args.out_csv}")
    print(summary[["quantity", "cpw_width_um", "cpw_gap_um", "n_lengths",
                   "intercept_frac_of_typical", "max_rel_fit_residual"]].to_string(index=False))

    bad = summary[summary["max_rel_fit_residual"] > 0.02]
    if len(bad):
        print("\nWARNING: >2% deviation from a straight line for:")
        print(bad[["quantity", "cpw_width_um", "cpw_gap_um", "max_rel_fit_residual"]]
             .to_string(index=False))
    else:
        print("\nAll points fit a straight line to within 2% -- "
              "C_total/length and L_total/integration_length are well-defined "
              "per-length quantities over this length range.")

    if not args.no_plots:
        make_plots(cap_df, ind_df, HERE)


if __name__ == "__main__":
    main()
