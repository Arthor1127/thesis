"""
Check mesh convergence of a run_amr.py run: read the per-AMR-iteration
outputs (iteration*/eig.csv + port-EPR.csv), identify the readout-resonator
mode in each, and track BOTH

  * the resonator frequency f_m           (target: < 0.5% between last two iters)
  * the coupling c derived from p         (target: < 5%   between last two iters)

The second check is the one EPR_C_EXTRACTION.md Sec 6.2 insists on: the
eigensolver's backward error says nothing about the field accuracy at the
3 um junction gaps, and c ~ sqrt(p).

Usage:  python check_convergence.py <run_dir> [--plot out.png]

<run_dir> holds epr_result.json; Palace output is at <run_dir>/*/outputFiles/.
"""
import argparse
import glob
import json
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import identify_resonator_mode, CJ_J_FF, JJ_PER_LOOP

PHI0_SQ_OVER_H_GHZ_NH = 163.4


def c_array_from_P(P_row, f_ghz, n_squids, L_J_nH):
    """Signed sum of per-loop c_k (same reformulation as extract_c_array.py)."""
    E_eff = PHI0_SQ_OVER_H_GHZ_NH / (2.0 * L_J_nH)
    pj = P_row
    tot = 0.0
    for k in range(n_squids):
        p = pj[2 * k] + pj[2 * k + 1]
        if p != 0:
            tot += math.copysign(0.5 * math.sqrt(abs(p) * f_ghz / (2.0 * E_eff)), p)
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--plot", default=None)
    args = ap.parse_args()

    meta = json.load(open(os.path.join(args.run_dir, "epr_result.json")))
    n_squids = meta["n_squids"]
    L_J_nH = meta["l_j_nh"]

    of = next(iter(glob.glob(os.path.join(args.run_dir, "*", "outputFiles"))
                   + glob.glob(os.path.join(args.run_dir, "outputFiles"))), None)
    if of is None:
        sys.exit(f"no outputFiles under {args.run_dir}")

    it_dirs = sorted(glob.glob(os.path.join(of, "iteration*")),
                     key=lambda d: int(d.split("iteration")[-1]))
    stages = it_dirs + [of]

    rows = []
    for d in stages:
        ep, ec = os.path.join(d, "port-EPR.csv"), os.path.join(d, "eig.csv")
        if not (os.path.exists(ep) and os.path.exists(ec)):
            continue
        P = pd.read_csv(ep).to_numpy()[:, 1:]
        fz = pd.read_csv(ec).to_numpy()[:, 1] * 1e9
        try:
            ri, _ = identify_resonator_mode(P, fz)
        except RuntimeError:
            continue
        f_ghz = float(fz[ri] / 1e9)
        sabs = float(np.abs(P[ri]).sum())
        c = c_array_from_P(P[ri], f_ghz, n_squids, L_J_nH) if P.shape[1] == JJ_PER_LOOP * n_squids else None
        rows.append({"stage": os.path.basename(d) or "final", "f_ghz": f_ghz,
                     "sum_abs_p": sabs, "c_array": c})

    # dedupe identical frequencies (final dir often duplicates the last iteration)
    uniq, seen = [], set()
    for r in rows:
        key = round(r["f_ghz"], 7)
        if key not in seen:
            seen.add(key)
            uniq.append(r)

    for i in range(1, len(uniq)):
        uniq[i]["df_pct"] = 100 * abs(uniq[i]["f_ghz"] - uniq[i - 1]["f_ghz"]) / uniq[i]["f_ghz"]
        if uniq[i]["c_array"] and uniq[i - 1]["c_array"]:
            uniq[i]["dc_pct"] = 100 * abs(uniq[i]["c_array"] - uniq[i - 1]["c_array"]) / abs(uniq[i]["c_array"])

    f_conv = len(uniq) >= 2 and uniq[-1].get("df_pct", 99) < 0.5
    c_conv = len(uniq) >= 2 and uniq[-1].get("dc_pct", 99) < 5.0

    print(f"\n{args.run_dir}  (L_J={L_J_nH}nH, n={n_squids})")
    print(f"{'stage':>12} {'f (GHz)':>11} {'df %':>8} {'sum|p_j|':>12} {'c_array':>12} {'dc %':>8}")
    for r in uniq:
        print(f"{r['stage']:>12} {r['f_ghz']:>11.5f} {r.get('df_pct', float('nan')):>8.3f} "
              f"{r['sum_abs_p']:>12.3e} "
              f"{(r['c_array'] if r['c_array'] is not None else float('nan')):>12.4e} "
              f"{r.get('dc_pct', float('nan')):>8.2f}")
    print(f"\nfrequency converged (<0.5%): {f_conv}")
    print(f"c_array   converged (<5%)  : {c_conv}")
    verdict = "CONVERGED" if (f_conv and c_conv) else "NOT converged - add AMR iterations / finer base mesh"
    print(f"VERDICT: {verdict}")

    if args.plot and len(uniq) >= 2:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        x = list(range(len(uniq)))
        fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
        ax[0].plot(x, [r["f_ghz"] for r in uniq], "o-")
        ax[0].set_xlabel("AMR stage"); ax[0].set_ylabel("resonator f (GHz)")
        ax[0].set_xticks(x); ax[0].set_xticklabels([r["stage"] for r in uniq], rotation=30)
        ax[0].set_title("frequency convergence")
        cvals = [r["c_array"] for r in uniq if r["c_array"] is not None]
        if cvals:
            ax[1].plot(x[:len(cvals)], cvals, "s-", color="C1")
            ax[1].set_xlabel("AMR stage"); ax[1].set_ylabel("c_array (rad)")
            ax[1].set_title("coupling convergence")
        plt.suptitle(f"{args.run_dir}  L_J={L_J_nH}nH n={n_squids}")
        plt.tight_layout(); plt.savefig(args.plot, dpi=110)
        print(f"plot -> {args.plot}")

    json.dump({"run_dir": args.run_dir, "l_j_nh": L_J_nH, "n_squids": n_squids,
               "stages": uniq, "frequency_converged": f_conv, "c_converged": c_conv},
              open(os.path.join(args.run_dir, "convergence.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
