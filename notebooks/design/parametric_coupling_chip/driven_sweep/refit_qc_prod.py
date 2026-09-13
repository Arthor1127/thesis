#!/usr/bin/env python
"""Batch refit for the production Qc sweep (out_qc_prod/) and the concluding
validation pair (out_final/).

run_qc_sim.py's shallow-dip guard (dip_depth_db > -3.0 -> skip fit) tripped
for EVERY point in both directories -- the dips are real but shallow
(-0.6 to -1.6dB net notch on a curving background), same pattern already
seen and fixed for cs14_cl300/cs10_cl400 earlier this session. Re-fits each
point around its own recorded dip_freq_ghz +-100MHz, the same recipe
validated in refit_convergence.py.

Flags (does NOT silently refit) any point whose recorded dip_freq_ghz sits
within WINDOW_HALF_GHZ of the sweep's own frequency edges (7.6/8.35GHz) --
for those the "dip" may just be the window boundary, not a real resonance,
and a +-100MHz refit window would run off the end of the data.

Run ON THE CLUSTER (never on the frontend):
    qsub submit_refit_qc_prod.sh
"""
import glob
import json
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo/data root, one level up from this subfolder
WINDOW_HALF_GHZ = 0.1
EDGE_MARGIN_GHZ = 0.05  # flag dips within this of fstart/fend as edge artifacts

# Sanity check added after the first production refit pass: run_qc_sim.py's
# own dip_depth_db is computed against the FULL 7.6-8.35GHz window, so a
# smoothly sloping background (unrelated to any resonance) reads as a
# multi-dB "dip" even when there is no real notch at dip_freq_ghz. Visual
# inspection of several npz traces (cs10_cl100, cs25_cl150) showed the
# fitted "resonance" was actually a ~0.003dB ripple riding on that slope --
# fit_transmission still returns a formally converged (Qr, Qc, fr), but it is
# fitting noise, not a physical notch, and gives absurd/wildly non-monotonic
# Qc vs coupling length as a result. LOCAL_DEPTH_DB_MIN filters this: only
# trust a refit whose dip is at least this deep relative to the LOCAL
# baseline (the +-WINDOW_HALF_GHZ window's own edges), not the global slope.
LOCAL_DEPTH_DB_MIN = 0.02


def refit_dir(results_dir, fit_transmission):
    rows = []
    edge_flagged = []
    for jf in sorted(glob.glob(os.path.join(results_dir, "*.json"))):
        with open(jf) as fh:
            r = json.load(fh)
        tag = r["tag"]
        npz_path = os.path.join(results_dir, f"{tag}_s21.npz")
        if not os.path.exists(npz_path):
            print(f"{tag}: no s21 npz found, skipping")
            continue

        f0 = r["dip_freq_ghz"]
        fstart, fend = r.get("fstart_ghz"), r.get("fend_ghz")
        near_edge = (fstart is not None and f0 - fstart < EDGE_MARGIN_GHZ) or \
                    (fend is not None and fend - f0 < EDGE_MARGIN_GHZ)
        if near_edge:
            print(f"{tag}: dip_freq_ghz={f0} is within {EDGE_MARGIN_GHZ}GHz of the "
                  f"sweep window edge [{fstart},{fend}] -- likely NOT the true "
                  f"resonance (window too narrow for this point). FLAGGED, not fit.")
            r["fr_ghz_refit"] = None
            r["Qr_refit"] = None
            r["Qc_refit"] = None
            r["Qi_refit"] = None
            r["edge_flagged"] = True
            r["local_depth_db"] = None
            r["weak_local_dip"] = None
            edge_flagged.append(tag)
            rows.append(r)
            continue

        d = np.load(npz_path)
        freqs = d["freq_hz"]
        S21 = d["S21"]
        f0_hz = f0 * 1e9
        mask = (freqs > f0_hz - WINDOW_HALF_GHZ * 1e9) & (freqs < f0_hz + WINDOW_HALF_GHZ * 1e9)
        if mask.sum() < 20:
            print(f"{tag}: not enough points in +-{WINDOW_HALF_GHZ}GHz window "
                  f"({mask.sum()}), skipping")
            continue

        mag_db = 20 * np.log10(np.abs(S21[mask]))
        n_edge = max(5, mask.sum() // 10)
        local_baseline = np.median(np.concatenate([mag_db[:n_edge], mag_db[-n_edge:]]))
        local_depth_db = local_baseline - mag_db.min()
        weak_local_dip = local_depth_db < LOCAL_DEPTH_DB_MIN
        if weak_local_dip:
            print(f"{tag}: local_depth_db={local_depth_db:.4f} < {LOCAL_DEPTH_DB_MIN} -- "
                  f"'dip' is background slope/noise, not a resolvable resonance. "
                  f"FLAGGED, not trusted (fit still computed for reference).")

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
        r["edge_flagged"] = False
        r["local_depth_db"] = local_depth_db
        r["weak_local_dip"] = weak_local_dip
        rows.append(r)
        tag_str = f"{tag} [WEAK]" if weak_local_dip else tag
        print(f"{tag_str}: fr={r['fr_ghz_refit']} GHz  Qr={Qr}  Qc={Qc}  Qi={Qi}  "
              f"local_depth_db={local_depth_db:.4f}")

    return rows, edge_flagged


def main():
    from qiskit_metal.analyses.em.transmission_fitting import fit_transmission

    all_rows = []
    for label, results_dir in [("out_qc_prod", os.path.join(HERE, "out_qc_prod", "results")),
                                ("out_final", os.path.join(HERE, "out_final", "results"))]:
        if not os.path.isdir(results_dir):
            print(f"{label}: no results dir at {results_dir}, skipping")
            continue
        print(f"\n=== refitting {label} ===")
        rows, edge_flagged = refit_dir(results_dir, fit_transmission)
        for r in rows:
            r["source"] = label
        all_rows.extend(rows)
        if edge_flagged:
            print(f"\n{label}: {len(edge_flagged)} point(s) flagged as window-edge "
                  f"artifacts (need a wider/shifted sweep window to re-run): {edge_flagged}")

    if not all_rows:
        raise SystemExit("nothing refit -- no results found")

    df = pd.DataFrame(all_rows)
    out_csv = os.path.join(HERE, "qc_prod_refit.csv")
    df.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")

    trusted = df[~df["edge_flagged"].fillna(False).astype(bool) &
                 ~df["weak_local_dip"].fillna(False).astype(bool)].copy()
    weak = df[df["weak_local_dip"].fillna(False).astype(bool)].copy()
    print(f"\n{len(trusted)}/{len(df)} points TRUSTED (real, locally-resolved notch, "
          f"local_depth_db >= {LOCAL_DEPTH_DB_MIN}dB).")
    print(f"{len(weak)} point(s) have a formally-converged fit but a LOCAL dip depth "
          f"< {LOCAL_DEPTH_DB_MIN}dB -- fitting background slope/noise, not trusted:")
    if len(weak):
        print(weak[["tag", "coupling_space_um", "coupling_length_um", "Qc_refit",
                    "local_depth_db"]].sort_values(
                    ["coupling_space_um", "coupling_length_um"]).to_string(index=False))
    print(f"{len(df) - len(trusted) - len(weak)} edge-flagged (window too narrow) or failed.")
    print("\nTRUSTED points:")
    print(trusted[["tag", "coupling_space_um", "coupling_length_um", "fr_ghz_refit",
              "Qc_refit", "Qi_refit", "local_depth_db"]].sort_values(
              ["coupling_space_um", "coupling_length_um"]).to_string(index=False))
    ok = trusted

    prod = ok[ok["source"] == "out_qc_prod"].copy()
    if len(prod) < 2:
        print("not enough production points for a Qc(cs, cl) plot, skipping figure")
        return

    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import thesis_style as ts
    ts.use()

    fig, ax = plt.subplots(figsize=ts.figsize(1))
    for cs, grp in prod.groupby("coupling_space_um"):
        grp = grp.sort_values("coupling_length_um")
        ax.semilogy(grp["coupling_length_um"], grp["Qc_refit"], "o-",
                    label=fr"$s={cs:g}\,\mu$m")
    ax.set_xlabel(r"coupling length $\ell$ ($\mu$m)")
    ax.set_ylabel(r"$Q_c$")
    ax.legend(fontsize=6, ncol=2)
    ts.finalize([ax])
    IMAGES_DIR = os.path.join(HERE, "images")
    os.makedirs(IMAGES_DIR, exist_ok=True)
    ts.save(fig, "qc_prod_refit", directory=IMAGES_DIR)
    print(f"\nwrote qc_prod_refit.{{png,pdf}} in {HERE}")


if __name__ == "__main__":
    main()
