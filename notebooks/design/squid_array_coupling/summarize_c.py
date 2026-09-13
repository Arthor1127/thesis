"""
Turn one or more run_amr.py epr_result.json files into the final c-extraction
summary: runs extract_c_array.py's analysis, prints a table, and (optionally)
plots c_array vs L_J, c_array vs n_squids, and the per-loop c_k profile.

Usage:
    python summarize_c.py campaign/*/epr_result.json conv_*/epr_result.json \
        --out c_summary.json --plotdir c_plots
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import extract_c_array as X


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results", nargs="+")
    ap.add_argument("--out", default="c_summary.json")
    ap.add_argument("--plotdir", default=None)
    args = ap.parse_args()

    runs = []
    for p in args.results:
        d = json.load(open(p))
        if d.get("status") != "success":
            print(f"SKIP {p}: {d.get('status')}")
            continue
        a = X.analyse_run(d)
        a["path"] = p
        a["amr_converged"] = d.get("amr", {}).get("frequency_converged_lt_0.5pct")
        runs.append(a)
    runs.sort(key=lambda r: (r["l_j_nh"], r["n_squids"]))

    print(f"\n{'L_J nH':>7} {'n':>4} {'f GHz':>9} {'sum|p_j|':>11} {'c_array':>12} "
          f"{'c_incoh':>12} {'dom loop':>9} {'conv':>6} {'elems':>10}")
    for r in runs:
        print(f"{r['l_j_nh']:>7.0f} {r['n_squids']:>4d} {r['f_ghz']:>9.4f} "
              f"{r['p_array_abs_sum']:>11.3e} {r['c_array']:>12.4e} "
              f"{r['c_array_incoherent']:>12.4e} {r['dominant_loop']['index']:>9d} "
              f"{str(r['amr_converged']):>6} {str(r['final_mesh_elements']):>10}")

    # L_J convergence (per n)
    by_n = {}
    for r in runs:
        by_n.setdefault(r["n_squids"], []).append(r)
    print("\n--- L_J convergence (EPR_C_EXTRACTION.md Sec 6.1 analog) ---")
    lj_out = []
    for n, g in sorted(by_n.items()):
        g = sorted(g, key=lambda r: r["l_j_nh"])
        if len(g) < 2:
            continue
        L = [r["l_j_nh"] for r in g]
        C = [abs(r["c_array"]) for r in g]
        spread = (max(C) - min(C)) / max(C)
        fit = X.fit_c_inf_L_loop(L, C)
        lj_out.append({"n_squids": n, "L_J_nH": L, "c_array": C, "spread": spread, "fit": fit})
        print(f"  n={n}: c_array {dict(zip(L, [round(c,6) for c in C]))}  spread={spread:.1%}"
              f"  ->  c_inf={fit['c_inf']:.4e}  L_loop={fit['L_loop_nH']*1000:.1f} pH"
              f"  ({'negligible, c~=c_inf' if spread<=0.05 else 'screening - use fit'})")

    # n_squids trend (per L_J)
    by_lj = {}
    for r in runs:
        by_lj.setdefault(r["l_j_nh"], []).append(r)
    print("\n--- coupling vs number of SQUIDs (fixed L_J) ---")
    nt_out = []
    for L_J, g in sorted(by_lj.items()):
        g = sorted(g, key=lambda r: r["n_squids"])
        if len(g) < 2:
            continue
        ns = [r["n_squids"] for r in g]
        cs = [r["c_array"] for r in g]
        mono = all(abs(cs[i]) <= abs(cs[i+1]) + 1e-12 for i in range(len(cs)-1))
        nt_out.append({"L_J_nH": L_J, "n_squids": ns, "c_array": cs, "increases_with_n": mono})
        print(f"  L_J={L_J}nH: " + "  ".join(f"n={n}:{c:.3e}" for n, c in zip(ns, cs))
              + f"   {'monotone up' if mono else 'NOT monotone'}")

    summary = {"runs": [{k: v for k, v in r.items() if k != "c_per_loop"} for r in runs],
               "lj_convergence": lj_out, "n_squids_trend": nt_out}
    json.dump(summary, open(args.out, "w"), indent=2)
    print(f"\nwrote {args.out}")

    if args.plotdir:
        os.makedirs(args.plotdir, exist_ok=True)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        # c vs L_J
        if lj_out:
            plt.figure(figsize=(6, 4.5))
            for e in lj_out:
                plt.plot(e["L_J_nH"], e["c_array"], "o-", label=f"n={e['n_squids']}")
                if e["fit"]:
                    plt.axhline(e["fit"]["c_inf"], ls=":", alpha=0.4)
            plt.xlabel("L_J (nH)"); plt.ylabel("c_array (rad)")
            plt.title("c vs junction inductance (dotted = c_inf fit)"); plt.legend()
            plt.tight_layout(); plt.savefig(f"{args.plotdir}/c_vs_LJ.png", dpi=120); plt.close()

        # c vs n_squids
        if nt_out:
            plt.figure(figsize=(6, 4.5))
            for e in nt_out:
                plt.plot(e["n_squids"], [abs(c) for c in e["c_array"]], "s-",
                         label=f"L_J={e['L_J_nH']}nH")
            plt.xlabel("number of SQUIDs"); plt.ylabel("|c_array| (rad)")
            plt.title("coupling vs array length"); plt.legend()
            plt.tight_layout(); plt.savefig(f"{args.plotdir}/c_vs_nsquids.png", dpi=120); plt.close()

        # per-loop c_k profile
        plt.figure(figsize=(8, 4.5))
        for r in runs:
            plt.plot(range(1, r["n_squids"] + 1), r["c_per_loop"], ".-",
                     label=f"L_J={r['l_j_nh']:.0f} n={r['n_squids']}")
        plt.xlabel("loop index (1 = far end, N = nearest the resonator)")
        plt.ylabel("c_k (rad)"); plt.title("per-loop coupling profile")
        plt.legend(fontsize=7); plt.tight_layout()
        plt.savefig(f"{args.plotdir}/c_per_loop.png", dpi=120); plt.close()
        print(f"plots -> {args.plotdir}/")


if __name__ == "__main__":
    main()
