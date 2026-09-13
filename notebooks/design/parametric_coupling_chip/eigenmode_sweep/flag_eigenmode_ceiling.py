#!/usr/bin/env python
"""Adds a reliability flag to qc_eigenmode.csv.

RETIRED the old Qc>82k/85k "likely_saturated"/"borderline" cliff cut. That cut
was bisected on the 1-port isolation test's loaded_Q (= fr/(2*Im f) off the
raw complex eigenfrequency), which does saturate at a ~4.5e7 noise floor past
Qc~85k. But the production Qc column here is NOT loaded_Q -- it's the kappa
route (fr / kappa_ext, kappa from port-Q.csv), and Phase C (job 1939013)
re-ran all 13 "stuck" grid points at solver_tol=1e-8 and found kappa-route Qc
reproducible to ~0.1-0.4% regardless of tolerance, with no ceiling up to the
highest value on the grid (~1e6). The cliff was a loaded_Q artifact
mis-applied to a different quantity. See HANDOFF.md / the plan for the full
trail. Every point is flagged "ok" here; independent cross-validation
(driven-S21, HFSS) still only covers Qc up to ~2e4 (see
qc_eigenmode_vs_driven.csv) -- that's a coverage gap, not a saturation flag.

Run ON THE CLUSTER: qsub submit_flag_ceiling.sh
"""
import os

import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo/data root, one level up from this subfolder


def main():
    path = os.path.join(HERE, "qc_eigenmode.csv")
    df = pd.read_csv(path)

    df["ceiling_flag"] = "ok"
    df.to_csv(path, index=False)

    print(f"Flagged {path} in place: all {len(df)} points -> ok "
          f"(loaded_Q cliff cut retired, see module docstring)")


if __name__ == "__main__":
    main()
