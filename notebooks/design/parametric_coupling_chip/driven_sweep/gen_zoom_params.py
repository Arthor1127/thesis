"""Build driven_sweep/params_qc_zoom.txt for the Phase-D narrow-window driven
sweep: one row `CS  CL  TAG  FR_SEED_GHZ` per stuck/borderline grid point, with
FR_SEED read from that point's eigenmode result (fr is trustworthy from the
eigenmode method even where Qc is not).

Run on the cluster (needs out_qc_eig/results/*.json).
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
EIG = os.path.join(REPO, "out_qc_eig", "results")

# The 13 points flagged likely_saturated / borderline in QC_STATUS.md.
STUCK = [
    (22.5, 150), (25, 300), (15, 150), (10, 100), (10, 150), (10, 250),
    (25, 250), (22.5, 200), (25, 150), (25, 100), (22.5, 250),
    (25, 350), (15, 250),
]

rows = []
for cs, cl in STUCK:
    tag = f"cs{cs:g}_cl{cl:g}"
    eig_json = os.path.join(EIG, f"{tag}_eig.json")
    if not os.path.exists(eig_json):
        print(f"WARNING: no eigenmode result for {tag} ({eig_json}) -- skipped")
        continue
    r = json.load(open(eig_json))
    fr = r.get("fr_ghz")
    if fr is None:
        print(f"WARNING: {tag} eigenmode result has no fr_ghz -- skipped")
        continue
    rows.append(f"{cs:g} {cl:g} {tag}_zoom {fr:.6f}")

out = os.path.join(HERE, "params_qc_zoom.txt")
with open(out, "w") as f:
    f.write("\n".join(rows) + "\n")
print(f"wrote {len(rows)} rows -> {out}")
print("\n".join(rows))
print(f"\nsubmit:  qsub -t 1-{len(rows)} driven_sweep/submit_qc_driven_zoom.sh")
