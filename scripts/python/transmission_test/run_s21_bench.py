"""
Self-contained driven (S21) run of the transmission_test design, for
benchmarking peak RAM + solve time on the cluster.

Mirrors build_and_driven.ipynb but calls the current build_design() (no args)
and writes NO ParaView fields (solns_to_save = 0).

Config via env vars: PALACE_BIN, NSLOTS, OUTDIR, FSTART, FEND, NPTS,
RES_MIN_UM, BG_MIN_UM, SOLVER_ORDER.
"""
import os
import sys
import time
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

OUTDIR     = os.environ.get("OUTDIR", "./s21_bench")
PALACE_BIN = os.environ["PALACE_BIN"]
NUM_CPUS   = int(os.environ.get("NSLOTS", "32"))
os.makedirs(OUTDIR, exist_ok=True)

FSTART_GHZ = float(os.environ.get("FSTART", "2.5"))
FEND_GHZ   = float(os.environ.get("FEND", "5.5"))
N_POINTS   = int(os.environ.get("NPTS", "6001"))
ADAPTIVE_TOL = 1e-3
ADAPTIVE_MAX_SAMPLES = 30

RES_MIN_UM   = float(os.environ.get("RES_MIN_UM", "3.0"))   # cpw_gap = 6 um -> ~half
BG_MIN_UM    = float(os.environ.get("BG_MIN_UM", "4.0"))
BG_MAX_UM    = 150.0
SOLVER_ORDER = int(os.environ.get("SOLVER_ORDER", "1"))

print(f"NUM_CPUS={NUM_CPUS}  sweep {FSTART_GHZ}-{FEND_GHZ} GHz x{N_POINTS}  "
      f"res_min={RES_MIN_UM}um  order={SOLVER_ORDER}", flush=True)

t0 = time.time()
from build_design import build_design

# qiskit-metal flips matplotlib to an interactive (Qt) backend on import;
# force it back so SQDMetal's headless mesh-preview PNG doesn't crash.
import matplotlib
matplotlib.use("Agg", force=True)

design = build_design()
design.rebuild()
print("components:", list(design.components.keys()), flush=True)
print(f"build+rebuild: {time.time() - t0:.1f} s", flush=True)

from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation

driven_options = {
    "dielectric_material": "silicon",
    "solns_to_save": 0,              # <-- no paraview / .pvd fields
    "solver_order": SOLVER_ORDER,
    "solver_tol": 1e-6,
    "solver_maxits": 500,
    "solver_initial_guess": True,
    "fillet_resolution": 12,
    "palace_dir": PALACE_BIN,
    "num_cpus": NUM_CPUS,
}

driven_sim = PALACE_Driven_Simulation(
    name="s21_bench", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", user_options=driven_options, create_files=True,
)
driven_sim.add_metallic(1)
driven_sim.add_ground_plane()
driven_sim.create_port_CPW_on_Launcher("LP1", 20e-6)
driven_sim.create_port_CPW_on_Launcher("LP2", 20e-6)
driven_sim.set_port_excitation(1)

driven_sim.fine_mesh_components(
    ["TL", "LP1", "LP2"], min_size=BG_MIN_UM * 1e-6, max_size=BG_MAX_UM * 1e-6,
    taper_dist_min=10e-6, metals_only=False)
driven_sim.fine_mesh_components(
    ["resonator1"], min_size=RES_MIN_UM * 1e-6, max_size=BG_MAX_UM * 1e-6,
    taper_dist_min=10e-6, metals_only=False)

fs, fe = FSTART_GHZ * 1e9, FEND_GHZ * 1e9
step = (fe - fs) / max(N_POINTS - 1, 1)
driven_sim.set_freq_values(freq_start=fs, freq_end=fe, freq_step=step)

t_prep = time.time()
driven_sim.prepare_simulation()
prep_s = time.time() - t_prep
print(f"MESH+PREP: {prep_s:.1f} s", flush=True)

# enable the reduced-order fast sweep
cfg_path = driven_sim._sim_config
with open(cfg_path) as f:
    cfg = json.load(f)
cfg["Solver"]["Driven"]["AdaptiveTol"] = ADAPTIVE_TOL
cfg["Solver"]["Driven"]["AdaptiveMaxSamples"] = ADAPTIVE_MAX_SAMPLES
with open(cfg_path, "w") as f:
    json.dump(cfg, f, indent=2)
print(f"PROM: AdaptiveTol={ADAPTIVE_TOL}, AdaptiveMaxSamples={ADAPTIVE_MAX_SAMPLES}", flush=True)

try:
    print("MESH SIZES:", driven_sim.retrieve_simulation_sizes(), flush=True)
except Exception as e:  # noqa
    print("retrieve_simulation_sizes failed:", e, flush=True)

t_solve = time.time()
driven_sim.run()
solve_s = time.time() - t_solve
print(f"SOLVE: {solve_s:.1f} s", flush=True)

s = driven_sim.retrieve_data()
freqs = np.asarray(s["freqs"])
S21 = np.asarray(s["S21"])
S11 = np.asarray(s.get("S11"))
np.savez(f"{OUTDIR}/s21.npz", freq_hz=freqs, S21=S21, S11=S11)
s21_db = 20 * np.log10(np.abs(S21))
i = int(np.argmin(s21_db))
print(f"|S21| min {s21_db[i]:.2f} dB @ {freqs[i] / 1e9:.5f} GHz   ({len(freqs)} pts reported)", flush=True)

print(f"\n=== TIMINGS ===  build {0:.0f}  mesh+prep {prep_s:.1f} s  solve {solve_s:.1f} s  "
      f"total {time.time() - t0:.1f} s", flush=True)
