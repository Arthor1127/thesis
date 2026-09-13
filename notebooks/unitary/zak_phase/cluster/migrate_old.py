#!/usr/bin/env python3
"""One-off: convert the first-run results/ (old schema) to results_chiC/ (new schema:
x_values / x_axis). No recompute."""
import glob
import json
import os

import numpy as np

SRC, DST = "results", "results_chiC"
for sd in sorted(glob.glob(os.path.join(SRC, "sweep_*"))):
    name = os.path.basename(sd)
    dst = os.path.join(DST, name)
    if os.path.exists(os.path.join(dst, "data.npz")):
        continue
    os.makedirs(dst, exist_ok=True)
    d = dict(np.load(os.path.join(sd, "data.npz")))
    d["x_values"] = d.pop("Chi_C")
    d["x_axis"] = "chi_C"
    np.savez(os.path.join(dst, "data.npz"), **d)
    p = json.load(open(os.path.join(sd, "params.json")))
    p["x_axis"] = "chi_C"
    p.setdefault("gamma_values", []); p.setdefault("omega_r_values", []); p.setdefault("chi_C_values", [])
    json.dump(p, open(os.path.join(dst, "params.json"), "w"), indent=2)
print(f"migrated {len(glob.glob(os.path.join(DST, 'sweep_*')))} sweeps -> {DST}")
