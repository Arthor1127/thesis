"""Eigenmode (lossy-port) extraction of Qc for one (coupling_space_um,
coupling_length_um) hanger resonator point -- an alternative to
run_qc_sim.py's driven-S21-notch-fit approach.

WHY: the production driven sweep (out_qc_prod/) showed every S21 notch in
this design's coupling regime is shallow (<0.25dB, several <0.02dB -- pure
background slope/noise), and the resulting Qc(coupling_length) was wildly
non-monotonic even after filtering out the worst noise-fit outliers. That is
inherent to fitting a Lorentzian on top of PALACE's absolute S21 error floor
for a WEAKLY coupled (= high Qc, small notch) resonator -- tighter solver
settings don't fix a signal-to-noise problem.

METHOD: give the SAME two CPW lumped ports used in the driven sim
(create_port_CPW_on_Launcher, default impedance_R=50 Ohm) to an EIGENMODE
simulation instead. A resistive (not open/PEC) port boundary makes the
eigenvalue problem lossy, so Palace returns a COMPLEX eigenfrequency
f = fr + i*fi for the resonator mode: fr is the resonance, and the standard
cavity-perturbation relation Q = fr/(2*fi) gives the mode's LOADED Q
directly from a single eigenmode solve -- no frequency sweep, no notch to
fit, no background curvature to fight.

Palace already separates this out per-port: per Palace's own docs, "By
specifying resistive lumped ports in the model, the mode coupling quality
factors will be computed as Qml = omega_m/kappa_ml", written to
port-Q.csv (kappa in Hz) alongside eig.csv's own overall loaded Q. Since our
ports are the ONLY loss channel besides silicon's tiny bulk loss tangent
(2.7e-6 -> Qi_bulk ~ 370000, negligible unless Qc approaches that scale),
kappa summed over both ports IS kappa_ext, and Qc = fr / kappa_ext.
SQDMetal's own PALACE_Eigenmode_Simulation.retrieve_mode_port_EPR() already
parses both eig.csv and port-Q.csv for exactly this (loaded_Q, kappa) pair
-- see Eigenmode_Simulation.py:446-495. No custom Palace config hacking
needed.

Reference: Ye, Wang, Liu (Tsinghua), "Electromagnetic Feature Extraction in
Superconducting Quantum Circuits: An Open-Source Finite-Element Workflow
Using Palace" (arXiv:2511.09041, Nov 2025) -- same Palace/GDSII workflow,
50-Ohm resistive lumped ports on the CPW feed, benchmarked resonator
frequencies (0.3%) and EXTERNAL COUPLING STRENGTHS (16%, 3/4 resonators)
against real cryogenic measurements. That is a materially stronger
precedent for trusting Palace's eigenmode-loss route for Qc than anything
found for the driven/S21 route (the SQDMetal paper itself never validated
driven S-parameter sims against anything).

MODE SELECTION: number_of_freqs modes are found starting from starting_freq
upward (Palace's shift-invert eigensolver). The physically real resonator
mode is the one with substantial kappa (real coupling to the lossy ports);
any spurious/artifact mode would have negligible kappa since it isn't
actually coupled to the ports. Picks argmax(kappa) as "primary".

Env: COUPLING_SPACE_UM, COUPLING_LENGTH_UM, PALACE_BIN, [NSLOTS=16],
     [OUTDIR=./out_qc_eig], [TAG], [STARTING_FREQ_GHZ=8.0], [NUM_FREQS=3],
     [SOLVER_ORDER=1], [SOLVER_TOL=1e-6], [SOLVER_MAXITS=100],
     [MESH_FACTOR=0.75] (matches the validated production driven-sim recipe)

Writes: {OUTDIR}/results/{TAG}.json

Cluster job only -- run via submit_qc_eigenmode.sh. Do not run this locally
(PALACE/FEM solves always go to the CNEA cluster).
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root, for readout_resonators etc.

COUPLING_SPACE_UM  = float(os.environ["COUPLING_SPACE_UM"])
COUPLING_LENGTH_UM = float(os.environ["COUPLING_LENGTH_UM"])
OUTDIR       = os.environ.get("OUTDIR", "./out_qc_eig")
PALACE_BIN   = os.environ["PALACE_BIN"]
NCPU         = int(os.environ.get("NSLOTS", "16"))
TAG          = os.environ.get("TAG", f"cs{COUPLING_SPACE_UM:g}_cl{COUPLING_LENGTH_UM:g}")
STARTING_FREQ_GHZ = float(os.environ.get("STARTING_FREQ_GHZ", "8.0"))
NUM_FREQS    = int(os.environ.get("NUM_FREQS", "3"))
SOLVER_ORDER = int(os.environ.get("SOLVER_ORDER", "1"))
SOLVER_TOL   = float(os.environ.get("SOLVER_TOL", "1e-8"))   # 1e-6 hit a hard Im(f)
                                                            # noise floor above Qc~85k
                                                            # (see ../QC_STATUS.md); 1e-8
                                                            # is also SQDMetal's own default
SOLVER_MAXITS = int(os.environ.get("SOLVER_MAXITS", "100"))
# See run_bareres_eigenmode_isolation.py: SQDMetal hardcodes
# Solver.Linear.Type="SuperLU"/KSPType="FGMRES"; set these to patch the config
# after prepare_simulation() (PalaceForCQED uses "Default"/"GMRES"). Empty = keep default.
LINEAR_TYPE  = os.environ.get("LINEAR_TYPE", "").strip()
KSP_TYPE     = os.environ.get("KSP_TYPE", "").strip()
MESH_FACTOR  = float(os.environ.get("MESH_FACTOR", "0.75"))  # matches validated driven-sim recipe
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

print(f"[{TAG}] coupling_space={COUPLING_SPACE_UM}um coupling_length={COUPLING_LENGTH_UM}um "
      f"ncpu={NCPU} starting_freq={STARTING_FREQ_GHZ}GHz num_freqs={NUM_FREQS}", flush=True)

t0 = time.time()
from qiskit_metal import designs  # noqa: F401  (forces the QtAgg->Agg reorder below)
import matplotlib
matplotlib.use("Agg", force=True)

from readout_resonators import build_design

design = build_design(COUPLING_SPACE_UM, COUPLING_LENGTH_UM)
print(f"geometry built ({time.time()-t0:.1f}s)  chip={design.chips.main.size.size_x} "
      f"x {design.chips.main.size.size_y}", flush=True)

from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation

eigen = PALACE_Eigenmode_Simulation(
    name=f"qceig_{TAG}", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "dielectric_material": "silicon",
        "starting_freq": STARTING_FREQ_GHZ * 1e9, "number_of_freqs": NUM_FREQS,
        "solns_to_save": NUM_FREQS,
        "solver_order": SOLVER_ORDER, "solver_tol": SOLVER_TOL, "solver_maxits": SOLVER_MAXITS,
        "palace_dir": PALACE_BIN, "num_cpus": NCPU,
    })
eigen.add_metallic(1)
eigen.add_ground_plane()
# SAME two ports as the driven sim (run_qc_sim.py), default impedance_R=50 Ohm
# -- this is what makes the eigenvalue problem lossy and gives Qc a meaning.
#
# PORT2_R_OHM: isolation knob added after the first real-design run
# (cs14_cl300_eig) found eig.csv's own loaded_Q (228) and its port-Q.csv-
# derived Qc (13562, summed over BOTH ports) disagreeing by ~60x -- while a
# bare single-lossy-port quarter-wave resonator (run_bareres_eigenmode_
# isolation.py) showed loaded_Q tracking the analytic Qc=(pi/4)(R/Z0) well
# (within 1-22%) on its own. That points at having TWO simultaneous lossy
# lumped ports (not the eigenmode method itself) as the likely culprit.
# Setting PORT2_R_OHM very large (near-open) turns LP2 back into a passive,
# effectively lossless boundary, isolating whether loaded_Q recovers its
# single-port reliability on the REAL design once only one lossy port remains.
PORT2_R_OHM = float(os.environ.get("PORT2_R_OHM", "50"))
eigen.create_port_CPW_on_Launcher("LP1", 30e-6)
eigen.create_port_CPW_on_Launcher("LP2", 30e-6, impedance_R=PORT2_R_OHM)

# Mesh: identical recipe to the validated production driven-sim setting
# (mesh_factor=0.75 -- see submit_qc_production.sh) so this is an
# apples-to-apples comparison, not a different-fidelity simulation.
lateral_min_um = 0.5 * MESH_FACTOR * min(20.0, 11.5, COUPLING_SPACE_UM)   # cpw_width, cpw_gap, coupling gap
res_min_um = max(lateral_min_um, 1.0)
bg_max_um = 200.0
eigen.fine_mesh_components(["TL1", "TL2", "LP1", "LP2"], min_size=8.0e-6, max_size=bg_max_um*1e-6,
                           taper_dist_min=16e-6, taper_dist_max=40e-6, metals_only=False)
eigen.fine_mesh_components(["readout_1", "CLT", "res_short"], min_size=res_min_um*1e-6,
                           max_size=bg_max_um*1e-6, taper_dist_min=30e-6, taper_dist_max=100e-6,
                           metals_only=False)
print(f"mesh: lateral_min={lateral_min_um:.3f}um res_min={res_min_um:.3f}um", flush=True)

tp = time.time()
eigen.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

if LINEAR_TYPE or KSP_TYPE:
    with open(eigen._sim_config) as f:
        cfg = json.load(f)
    if LINEAR_TYPE:
        cfg["Solver"]["Linear"]["Type"] = LINEAR_TYPE
    if KSP_TYPE:
        cfg["Solver"]["Linear"]["KSPType"] = KSP_TYPE
    with open(eigen._sim_config, "w") as f:
        json.dump(cfg, f, indent=2)
    print(f"PATCHED Solver.Linear -> {cfg['Solver']['Linear']}", flush=True)

ts = time.time()
eigen.run()
solve_s = time.time() - ts
print(f"SOLVE: {solve_s:.1f} s", flush=True)

eig_csv_path = eigen._output_data_dir + "/eig.csv"
if not os.path.exists(eig_csv_path):
    raise SystemExit(f"FAILED: run() completed without raising, but {eig_csv_path} was never "
                     "written (MPI job likely killed mid-run -- check for OOM/signal 9 above).")

port_epr = eigen.retrieve_mode_port_EPR()
mesh_sizes = eigen.retrieve_simulation_sizes()

freqs_hz = port_epr["eigenfrequencies"]
loaded_Q = port_epr["loaded_Q"]
# BUG WORKAROUND: SQDMetal's retrieve_mode_port_EPR_from_file() docstring
# claims 'kappa' is returned in Hz, but Palace's own port-Q.csv header is
# explicit: "kappa_ext[n] (GHz)" -- SQDMetal reads that column straight
# through with no unit conversion, so the array it hands back is actually in
# GHz. Confirmed directly against port-Q.csv from the cs14_cl300_eig run
# (2026-09-08): un-corrected this produced a nonsense Qc ~1e13. x1e9 here.
kappa_hz = port_epr["kappa"] * 1e9

import numpy as np
n_modes = len(freqs_hz)
if n_modes == 0:
    raise SystemExit(f"FAILED: zero converged modes (starting_freq={STARTING_FREQ_GHZ}GHz, "
                     f"number_of_freqs={NUM_FREQS}).")

modes = []
for i in range(n_modes):
    k = float(kappa_hz[i]) if np.isfinite(kappa_hz[i]) else None
    f_ghz = float(freqs_hz[i]) / 1e9
    Qc_i = f_ghz * 1e9 / k if (k is not None and k > 0) else None
    modes.append({"m": i + 1, "f_ghz": f_ghz, "loaded_Q": float(loaded_Q[i]),
                 "kappa_hz": k, "Qc": Qc_i})
    print(f"  mode {i+1}: f={f_ghz:.5f}GHz loaded_Q={loaded_Q[i]:.4e} "
         f"kappa={k}Hz Qc={Qc_i}", flush=True)
    # SANITY CHECK: if the port resistors are the only real loss channel
    # (silicon's bulk loss tangent alone gives Qi~370000, negligible unless
    # Qc approaches that), eig.csv's own loaded_Q should equal Qc_i to
    # within the small dielectric correction. First run (cs14_cl300, mode
    # 1) found loaded_Q=228 vs Qc_i(from port-Q.csv)=13562 -- a ~60x
    # internal Palace inconsistency, NOT explained by dielectric loss.
    # Flag it loudly rather than silently trust either number.
    if Qc_i is not None and loaded_Q[i] > 0:
        ratio = loaded_Q[i] / Qc_i
        if not (0.5 < ratio < 2.0):
            print(f"    WARNING: mode {i+1} loaded_Q/Qc_i = {ratio:.4f}, far from 1 -- "
                 f"eig.csv's own loaded_Q and port-Q.csv-derived Qc disagree. This is a "
                 f"real internal Palace inconsistency seen before, not just a units bug "
                 f"-- do not trust either number for this mode without investigating.",
                 flush=True)

# Primary/resonator mode = closest to STARTING_FREQ_GHZ.
# WAS: argmax(kappa) ("a spurious mode would have negligible kappa") -- WRONG,
# confirmed on cs10_cl400_eig/cs20_cl200_eig: a ~11.37GHz mode (clearly NOT
# our ~8.2GHz resonator -- likely a feedline/box mode) had a LARGER kappa
# than the real resonator mode in both cases, so argmax(kappa) silently
# reported the wrong mode's fr/Qc. Proximity to starting_freq is a much more
# direct, physically-grounded selector: we already know roughly where our
# resonator sits (from the driven fits), Palace's shift-invert eigensolver
# searches outward from starting_freq, so the intended mode should simply be
# the closest one to it.
primary = min(modes, key=lambda m: abs(m["f_ghz"] - STARTING_FREQ_GHZ))

row = {
    "tag": TAG, "coupling_space_um": COUPLING_SPACE_UM, "coupling_length_um": COUPLING_LENGTH_UM,
    "port2_r_ohm": PORT2_R_OHM,
    "fr_ghz": primary["f_ghz"] if primary else None,
    "Qc": primary["Qc"] if primary else None,
    "loaded_Q": primary["loaded_Q"] if primary else None,
    "kappa_hz": primary["kappa_hz"] if primary else None,
    "modes": modes,
    "starting_freq_ghz": STARTING_FREQ_GHZ, "number_of_freqs": NUM_FREQS, "num_cpus": NCPU,
    "solver_order": SOLVER_ORDER, "solver_tol": SOLVER_TOL, "solver_maxits": SOLVER_MAXITS,
    "linear_type": LINEAR_TYPE or "SuperLU", "ksp_type": KSP_TYPE or "FGMRES",
    "mesh_factor": MESH_FACTOR,
    "mesh_elements": mesh_sizes.get("MeshElements"),
    "solve_s": solve_s, "wall_s": time.time() - t0,
}
with open(f"{OUTDIR}/results/{TAG}.json", "w") as f:
    json.dump(row, f, indent=2)

print(json.dumps(row, indent=2), flush=True)

import glob
for pat in (f"{eigen._output_data_dir}/../*.msh", f"{eigen._output_data_dir}/paraview"):
    for p in glob.glob(pat):
        os.system(f"rm -rf '{p}'")
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
