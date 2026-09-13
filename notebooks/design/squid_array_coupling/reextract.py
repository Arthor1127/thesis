"""
Re-derive epr_result.json from a completed run's Palace output, WITHOUT
re-solving. Reads outputFiles/port-EPR.csv + eig.csv, re-identifies the
readout-resonator mode with build_design.identify_resonator_mode (sum|p_j|),
and rewrites <run_dir>/epr_result.json with the resonator mode's full
per-junction participation vector.

Usage:
    python reextract.py <run_dir> [L_J_nH]

<run_dir> is the directory that contains epr_result.json; the Palace output
is looked for at <run_dir>/*/outputFiles/ (or <run_dir>/outputFiles/).
"""
import glob
import json
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import (identify_resonator_mode, N_SQUIDS, L_J_NH, CJ_J_FF,
                          RES_LENGTH_MM, JJ_PER_LOOP)

run_dir = sys.argv[1]
L_J_nH = float(sys.argv[2]) if len(sys.argv) > 2 else L_J_NH

cands = (glob.glob(os.path.join(run_dir, "outputFiles"))
         + glob.glob(os.path.join(run_dir, "*", "outputFiles")))
out = next((d for d in cands if os.path.exists(os.path.join(d, "port-EPR.csv"))), None)
if out is None:
    sys.exit(f"no outputFiles/port-EPR.csv found under {run_dir}")
print(f"reading {out}")

P = pd.read_csv(os.path.join(out, "port-EPR.csv")).to_numpy()[:, 1:]   # (modes, ports)
eig = pd.read_csv(os.path.join(out, "eig.csv")).to_numpy()
freqs_hz = eig[:, 1] * 1e9
Q = eig[:, 3]
num_modes, n_ports = P.shape

p_signed = P.sum(axis=1)
p_abs = np.abs(P).sum(axis=1)
res_idx, res_info = identify_resonator_mode(P, freqs_hz)

omega0 = 2.0 * math.sqrt(1.0 / (L_J_nH * 1e-9 * CJ_J_FF * 1e-15))

resonator = {
    "m": int(res_idx + 1),
    "f_ghz": float(freqs_hz[res_idx] / 1e9),
    "Q": float(Q[res_idx]),
    "p_array_signed": float(p_signed[res_idx]),
    "p_array_abs_sum": float(p_abs[res_idx]),
    "p_per_junction": [float(x) for x in P[res_idx]],
}
if n_ports == JJ_PER_LOOP * N_SQUIDS:
    resonator["p_per_loop"] = [
        float(P[res_idx, 2 * i] + P[res_idx, 2 * i + 1]) for i in range(N_SQUIDS)
    ]

result = {
    "status": "success",
    "reextracted": True,
    "n_squids": N_SQUIDS,
    "l_j_nh": L_J_nH,
    "cj_j_ff": CJ_J_FF,
    "res_length_mm": RES_LENGTH_MM,
    "n_inductive_ports": int(n_ports),
    "resonator_id": res_info,
    "all_modes": [{
        "m": i + 1,
        "f_ghz": float(freqs_hz[i] / 1e9),
        "Q": float(Q[i]),
        "p_array_signed": float(p_signed[i]),
        "p_array_abs_sum": float(p_abs[i]),
        "kind": "resonator" if i == res_idx else "phonon/box",
    } for i in range(num_modes)],
    "resonator_mode": resonator,
    "phonon_dispersion": {
        "formula": "omega_j = 2*sqrt(1/(L_J*C)) * sin(j*pi/(2*(N+1)))  [acoustic limit, C = C_g]",
        "omega0_over_2pi_ghz_using_CJ": omega0 / (2 * math.pi) / 1e9,
        "observed_phonon_ghz": [
            float(freqs_hz[i] / 1e9) for i in range(num_modes) if i != res_idx
        ],
    },
}
dst = os.path.join(run_dir, "epr_result.json")
json.dump(result, open(dst, "w"), indent=2)

print(f"\nmodes: " + ", ".join(
    f"{freqs_hz[i]/1e9:.3f}GHz(|p|={p_abs[i]:.1e})" for i in range(num_modes)))
print(f"RESONATOR: m={res_idx+1}, f={freqs_hz[res_idx]/1e9:.4f} GHz, "
      f"sum|p_j|={p_abs[res_idx]:.2e}")
print(f"wrote {dst}")
