#!/usr/bin/env python
"""run_qc_sim.py's shallow-dip guard (dip_depth_db > -3.0 -> skip fit) tripped
for every point in the convergence sweep, same as it did earlier this session
for cs14_cl300 at the wide 7.6-8.35GHz window: a real but background-
dominated shallow dip (net ~0.3dB notch on a curving baseline), which needs a
narrowed window around the dip before fit_transmission can isolate it from
the background. Re-fits every point in out_conv/results/ using each point's
own recorded dip_freq_ghz +-100MHz, exactly the manual recipe already
validated for this same design point earlier in this session.

Run ON THE CLUSTER (same convention as aggregate_convergence.py):
    python refit_convergence.py
"""
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo/data root, one level up from this subfolder
WINDOW_HALF_GHZ = 0.1


def main():
    from qiskit_metal.analyses.em.transmission_fitting import fit_transmission

    rows = []
    for jf in sorted(glob.glob(os.path.join(HERE, "out_conv", "results", "*.json"))):
        with open(jf) as fh:
            r = json.load(fh)
        tag = r["tag"]
        npz_path = os.path.join(HERE, "out_conv", "results", f"{tag}_s21.npz")
        if not os.path.exists(npz_path):
            print(f"{tag}: no s21 npz found, skipping")
            continue
        d = np.load(npz_path)
        freqs = d["freq_hz"]; S21 = d["S21"]
        f0 = r["dip_freq_ghz"] * 1e9
        mask = (freqs > f0 - WINDOW_HALF_GHZ * 1e9) & (freqs < f0 + WINDOW_HALF_GHZ * 1e9)
        if mask.sum() < 20:
            print(f"{tag}: not enough points in +-{WINDOW_HALF_GHZ}GHz window "
                 f"({mask.sum()}), skipping")
            continue
        try:
            fit, _ = fit_transmission(freqs[mask], S21[mask], plot=False)
            Qr, Qc, fr = float(fit["Qr"]), float(fit["Qc"]), float(fit["fr"])
            Qi = 1.0 / (1.0 / Qr - 1.0 / Qc) if Qr != Qc else float("inf")
        except Exception as e:
            print(f"{tag}: refit FAILED -- {type(e).__name__}: {e}")
            Qr = Qc = Qi = fr = None

        r["fr_ghz_refit"] = fr / 1e9 if fr is not None else None
        r["Qr_refit"] = Qr
        r["Qc_refit"] = Qc
        r["Qi_refit"] = Qi
        rows.append(r)
        print(f"{tag}: fr={r['fr_ghz_refit']} GHz  Qr={Qr}  Qc={Qc}  Qi={Qi}")

    if not rows:
        raise SystemExit("nothing refit -- no results found")

    df = pd.DataFrame(rows)
    order = ["baseline", "tol_tighter", "tol_1e-5", "tol_1e-4", "tol_1e-4_maxits50",
            "tol_1e-3", "mesh_finer", "mesh_1.5x", "mesh_2x", "mesh_3x",
            "ams_more", "ams_200", "ams_120", "ams_80"]
    df["tag"] = pd.Categorical(df["tag"], categories=[t for t in order if t in df["tag"].values],
                               ordered=True)
    df = df.sort_values("tag")
    out_csv = os.path.join(HERE, "convergence_refit.csv")
    df.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")

    if "baseline" in df["tag"].values:
        base = df[df["tag"] == "baseline"].iloc[0]
        print(f"\nbaseline (refit): fr={base['fr_ghz_refit']:.5f} GHz, Qc={base['Qc_refit']:.1f}, "
             f"Qr={base['Qr_refit']:.1f}")
        df["fr_dev_pct"] = 100 * (df["fr_ghz_refit"] - base["fr_ghz_refit"]) / base["fr_ghz_refit"]
        df["Qc_dev_pct"] = 100 * (df["Qc_refit"] - base["Qc_refit"]) / base["Qc_refit"]
        print(df[["tag", "solver_tol", "solver_maxits", "mesh_factor", "adaptive_max_samples",
                  "fr_ghz_refit", "Qc_refit", "fr_dev_pct", "Qc_dev_pct",
                  "solve_s"]].to_string(index=False))
        flagged = df[(df["fr_dev_pct"].abs() > 1.0) | (df["Qc_dev_pct"].abs() > 15.0)]
        if len(flagged):
            print(f"\nWARNING: {len(flagged)} point(s) deviate from baseline by >1% (fr) or "
                 f">15% (Qc):")
            print(flagged[["tag", "fr_dev_pct", "Qc_dev_pct"]].to_string(index=False))
        else:
            print("\nAll refit points agree with baseline within 1% (fr) / 15% (Qc).")

    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import thesis_style as ts
    ts.use()

    fig, axs = plt.subplots(2, 1, figsize=ts.figsize(1, ratio=1.2))
    x = np.arange(len(df))
    labels = df["tag"].astype(str).tolist()
    base_fr = df[df["tag"] == "baseline"]["fr_ghz_refit"].iloc[0] if "baseline" in df["tag"].values else None
    base_Qc = df[df["tag"] == "baseline"]["Qc_refit"].iloc[0] if "baseline" in df["tag"].values else None

    if base_fr is not None:
        axs[0].axhline(base_fr, ls="--", lw=0.8, color="gray")
    axs[0].plot(x, df["fr_ghz_refit"], "o", color=ts.CYCLE[0])
    axs[0].set_ylabel(r"$f_r$ (GHz, refit)")

    if base_Qc is not None:
        axs[1].axhline(base_Qc, ls="--", lw=0.8, color="gray")
    axs[1].semilogy(x, df["Qc_refit"], "o", color=ts.CYCLE[1])
    axs[1].set_ylabel(r"$Q_c$ (refit)")
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(labels, rotation=60, ha="right", fontsize=6)

    ts.finalize(axs)
    ts.panel_labels(axs)
    IMAGES_DIR = os.path.join(HERE, "images")
    os.makedirs(IMAGES_DIR, exist_ok=True)
    ts.save(fig, "convergence_refit", directory=IMAGES_DIR)
    print(f"\nwrote convergence_refit.{{png,pdf}} in {HERE}")


if __name__ == "__main__":
    main()
