#!/usr/bin/env python
"""One-time fixup: 3 result files in out_qc_eig/results/ were written before
run_qc_eigenmode.py's mode-selection bug was fixed (argmax(kappa) -> closest
to starting_freq_ghz) and, for cs14_cl300_eig specifically, also before the
port-Q.csv units bug was fixed (GHz mislabeled as Hz). No re-simulation
needed -- each file's own "modes" list already has every mode's data; this
just re-derives the top-level fr_ghz/loaded_Q/kappa_hz/Qc fields using the
CURRENT (correct) selection rule, in place.

Run ON THE CLUSTER (postprocessing only, but still not the frontend):
    qsub submit_fixup_stale.sh
"""
import glob
import json
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo/data root, one level up from this subfolder
STALE_TAGS = {"cs14_cl300_eig", "cs14_cl300_1port", "cs20_cl200_eig"}


def main():
    results_dir = os.path.join(HERE, "out_qc_eig", "results")
    fixed = []
    for jf in sorted(glob.glob(os.path.join(results_dir, "*.json"))):
        with open(jf) as f:
            r = json.load(f)
        if r.get("tag") not in STALE_TAGS:
            continue
        modes = r.get("modes", [])
        if not modes:
            print(f"{r['tag']}: no 'modes' list, skipping")
            continue
        starting_freq = r.get("starting_freq_ghz", 8.0)

        # Re-derive kappa_hz/Qc for every mode too, in case this file
        # predates the units fix (GHz mislabeled as Hz -- symptom: Qc >> 1e6).
        for m in modes:
            if m.get("Qc") and m["Qc"] > 1e8 and m.get("kappa_hz"):
                m["kappa_hz"] = m["kappa_hz"] * 1e9
                m["Qc"] = m["f_ghz"] * 1e9 / m["kappa_hz"]

        primary = min(modes, key=lambda m: abs(m["f_ghz"] - starting_freq))
        old = (r.get("fr_ghz"), r.get("Qc"))
        r["fr_ghz"] = primary["f_ghz"]
        r["loaded_Q"] = primary["loaded_Q"]
        r["kappa_hz"] = primary["kappa_hz"]
        r["Qc"] = primary["Qc"]
        r["modes"] = modes
        r["_fixed_up"] = True

        with open(jf, "w") as f:
            json.dump(r, f, indent=2)
        print(f"{r['tag']}: fr/Qc {old} -> ({r['fr_ghz']}, {r['Qc']})")
        fixed.append(r["tag"])

    print(f"\nFixed {len(fixed)}/{len(STALE_TAGS)} expected stale files: {fixed}")
    missing = STALE_TAGS - set(fixed)
    if missing:
        print(f"WARNING: expected but not found/fixed: {missing}")


if __name__ == "__main__":
    main()
