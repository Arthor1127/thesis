"""Driven S21 simulation of one (coupling_space_um, coupling_length_um) hanger
resonator point (see readout_resonators.py), fitted for resonance frequency
and quality factors with qiskit_metal's own phi-RM fitter
(qiskit_metal.analyses.em.transmission_fitting.fit_transmission) --
the same PALACE_Driven_Simulation + lumped-CPW-port pattern as
scripts/python/transmission_test/run_s21_qw.py.

Env: COUPLING_SPACE_UM, COUPLING_LENGTH_UM, PALACE_BIN, [NSLOTS=16],
     [OUTDIR=./out_qc], [TAG], [FSTART_GHZ=7.6], [FEND_GHZ=8.35], [NPTS=12001],
     [ADAPTIVE_TOL=1e-3], [ADAPTIVE_MAX_SAMPLES=300], [SOLVER_ORDER=1]

NOTE (this session): readout_resonators.py's snake geometry was rebuilt from
scratch since this window was first tuned (bigger fillets, single-fold
routing, COUPLING_BUDGET_UM 550->700, a length-safety floor on leg lengths)
-- the routed resonator length changed as a side effect. An eigenmode check
with NO ports (build_and_eigenmode.py, run across 4 different feedline
lengths for robustness) found the bare resonator mode sitting at a very
consistent 8.25-8.28GHz (30MHz spread), not the original 8.0GHz target. The
window below is re-centered on that -- verify it still finds a real,
non-shallow dip before trusting any point's fit.
Writes: {OUTDIR}/results/{TAG}.json, {OUTDIR}/results/{TAG}_s21.npz,
        {OUTDIR}/results/{TAG}_fit.png

Cluster job only -- run via submit_qc_sweep.sh. Do not run this locally
(PALACE/FEM solves always go to the CNEA cluster).

Mesh sizing follows this repo's "half the smallest lateral feature" rule
(palace-mesh-resolution-and-amr): the smallest feature that actually
determines Qc is the coupling_space_um gap itself, so the coupling region's
fine-mesh target is tied to it, not to a fixed constant.

FREQUENCY WINDOW -- this is the one setting most likely to silently ruin a
sweep point, so read this before changing it. The first pilot (FSTART=7.0,
FEND=9.0, ADAPTIVE_MAX_SAMPLES=40) came back with a nonsense fit (Qr=Qi=32 --
for a lossless PEC/vacuum sim, Qi should be enormous, not comparable to Qr).
Plotting the raw |S21|/phase showed why: the curve was essentially FLAT
across the whole 2GHz window (a ~0.07dB wiggle, no phase wrap) -- PROM's
adaptive sampler, spread over 2GHz with only 40 samples (~50MHz average
spacing), never happened to land near the actual (likely much narrower)
resonance, so its greedy refinement had no signal to zoom in on and converged
to an uninformative smooth model. `fit_transmission` then dutifully fit a
Lorentzian to noise.

Confirmed by checking PALACE's own config handling
(SQDMetal.PALACE.Frequency_Driven_Simulation.set_freq_values writes
MinFreq/MaxFreq/FreqStep): FreqStep/NPTS is ONLY the output grid PROM's cheap
reduced-order model gets evaluated/written on -- it has NO effect on which
frequencies get real (expensive) HDM training solves; that is chosen
adaptively within [MinFreq, MaxFreq] up to AdaptiveMaxSamples alone. So a
finer step over the SAME wide window (e.g. 0.2 MHz instead of 0.5 MHz) cannot
fix a missed resonance -- only shrinking the window (concentrating the same
sample budget where the resonance actually lives) or raising
AdaptiveMaxSamples can. The window here is deliberately narrow AND biased
below the bare design target of 8.0 GHz: capacitive loading from coupling to
the feedline pulls f_r down somewhat (a real, expected effect for a hanger
resonator). The 120 MHz window (7.90-8.02 GHz) still missed it for the
tightest-coupling pilot point (coupling_space=10um, coupling_length=100um --
the strongest pull in the whole grid): |S21| was monotonic and still moving
at the 7.90 GHz edge, meaning the true dip sits below it. Widened to
FSTART/FEND=7.5/8.02 GHz (520 MHz) with ADAPTIVE_MAX_SAMPLES=300 (~1.7 MHz
average training-sample spacing) for the next attempt. Still a
physically-motivated guess, not a guarantee -- ALWAYS check dip_depth_db in
the result (and/or the saved *_fit.png) before trusting Qr/Qc/Qi for any
point: a shallow dip (order -1dB, as the failed pilots showed) means the
window missed the true resonance for that (coupling_space_um,
coupling_length_um) point and it needs a re-run with a shifted/widened
window, not a number to plot -- run_qc_sim.py itself now skips fit_transmission
entirely (Qr/Qc/Qi=None) when dip_depth_db > -3dB, since fitting a Lorentzian
to a flat/monotonic trace previously crashed curve_fit outright (ValueError:
initial guess outside provided bounds) after burning a full ~45-minute solve.
If several grid points come back shallow in the same direction (all above or
all below the window), that tells you which way to shift it.

NPTS=12001 (0.01 MHz step over the 120 MHz window) is a deliberately fine
OUTPUT grid for this diagnostic pilot, to see the lineshape clearly once
AdaptiveMaxSamples actually finds it -- this raises "Postprocessing" cost
(see the failed pilot's timing breakdown: 664s of 844s solve time was
Postprocessing, scaling roughly linearly with point count) well above the
NPTS=500-3000 range that's fine for routine sweep points once the window is
validated; drop it back down for the full 35-point sweep.
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root, for readout_resonators etc.

COUPLING_SPACE_UM  = float(os.environ["COUPLING_SPACE_UM"])
COUPLING_LENGTH_UM = float(os.environ["COUPLING_LENGTH_UM"])
OUTDIR       = os.environ.get("OUTDIR", "./out_qc")
PALACE_BIN   = os.environ["PALACE_BIN"]
NCPU         = int(os.environ.get("NSLOTS", "16"))
TAG          = os.environ.get("TAG", f"cs{COUPLING_SPACE_UM:g}_cl{COUPLING_LENGTH_UM:g}")
FSTART_GHZ   = float(os.environ.get("FSTART_GHZ", "7.6"))
FEND_GHZ     = float(os.environ.get("FEND_GHZ", "8.35"))
N_POINTS     = int(os.environ.get("NPTS", "12001"))
ADAPTIVE_TOL = float(os.environ.get("ADAPTIVE_TOL", "1e-3"))
ADAPTIVE_MAX_SAMPLES = int(os.environ.get("ADAPTIVE_MAX_SAMPLES", "300"))
SOLVER_ORDER = int(os.environ.get("SOLVER_ORDER", "1"))
SOLVER_TOL = float(os.environ.get("SOLVER_TOL", "1e-6"))
SOLVER_MAXITS = int(os.environ.get("SOLVER_MAXITS", "500"))
MESH_FACTOR = float(os.environ.get("MESH_FACTOR", "1.0"))  # multiplies the coupling-region min_size
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

print(f"[{TAG}] coupling_space={COUPLING_SPACE_UM}um coupling_length={COUPLING_LENGTH_UM}um "
      f"ncpu={NCPU} sweep {FSTART_GHZ}-{FEND_GHZ}GHz x{N_POINTS}", flush=True)

t0 = time.time()
from qiskit_metal import designs  # noqa: F401  (forces the QtAgg->Agg reorder below)
import matplotlib
matplotlib.use("Agg", force=True)

from readout_resonators import build_design

design = build_design(COUPLING_SPACE_UM, COUPLING_LENGTH_UM)
print(f"geometry built ({time.time()-t0:.1f}s)  chip={design.chips.main.size.size_x} "
      f"x {design.chips.main.size.size_y}", flush=True)

from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation

driven = PALACE_Driven_Simulation(
    name=f"qc_{TAG}", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "dielectric_material": "silicon", "solns_to_save": 0,
        "solver_order": SOLVER_ORDER, "solver_tol": SOLVER_TOL, "solver_maxits": SOLVER_MAXITS,
        "solver_initial_guess": True, "fillet_resolution": 12,
        "palace_dir": PALACE_BIN, "num_cpus": NCPU,
    })
driven.add_metallic(1)
driven.add_ground_plane()
driven.create_port_CPW_on_Launcher("LP1", 30e-6)
driven.create_port_CPW_on_Launcher("LP2", 30e-6)
driven.set_port_excitation(1)

# base mesh: half the smallest lateral feature (palace-mesh-resolution-and-amr) --
# the coupling gap itself is what determines Qc, so it sets the fine-mesh floor.
AMR_ITERS = int(os.environ.get("AMR_ITERS", "0"))
if AMR_ITERS > 0 and MESH_FACTOR != 1.0:
    print(f"WARNING: AMR_ITERS>0 but MESH_FACTOR={MESH_FACTOR} != 1.0 -- forcing "
         f"MESH_FACTOR=1.0 for the base mesh. AMR refines from an existing mesh's own "
         f"error indicator; it cannot rescue a base mesh too coarse to resolve the "
         f"coupling gap to begin with (coarsening it 1.5x alone already gave a "
         f"qualitatively wrong Qc, +208%, in this session's manual convergence sweep).",
         flush=True)
    MESH_FACTOR = 1.0
lateral_min_um = 0.5 * MESH_FACTOR * min(20.0, 11.5, COUPLING_SPACE_UM)   # cpw_width, cpw_gap, coupling gap
res_min_um = max(lateral_min_um, 1.0)
bg_max_um = 200.0
driven.fine_mesh_components(["TL1", "TL2", "LP1", "LP2"], min_size=8.0e-6, max_size=bg_max_um*1e-6,
                            taper_dist_min=16e-6, taper_dist_max=40e-6, metals_only=False)
driven.fine_mesh_components(["readout_1", "CLT", "res_short"], min_size=res_min_um*1e-6,
                            max_size=bg_max_um*1e-6, taper_dist_min=30e-6, taper_dist_max=100e-6,
                            metals_only=False)
print(f"mesh: lateral_min={lateral_min_um:.3f}um res_min={res_min_um:.3f}um", flush=True)

# Solution-based AMR (spatial mesh refinement during the solve), NOT the same
# thing as PROM's own frequency-domain adaptive sampling above. AMR refines
# based on an error indicator computed on the EXISTING mesh, so it can only
# ever sharpen a starting point that already resolves the geometry reasonably
# -- it cannot rescue a base mesh too coarse to represent the coupling gap's
# field at all. AMR_ITERS was already read above (before the base mesh was
# built) so MESH_FACTOR could be force-corrected to 1.0 first if needed.
AMR_TOL = float(os.environ.get("AMR_TOL", "1e-2"))
# Default capped (NOT 0/unbounded): mesh_finer alone (a uniform, one-shot 2x
# refinement of this same base mesh) needed 103.6GB -- AMR refining
# repeatedly and specifically at the gap where Qc is most sensitive could
# plausibly run away similarly or worse if left unbounded.
AMR_MAX_DOFS = int(os.environ.get("AMR_MAX_DOFS", "1500000"))
AMR_DORFLER = float(os.environ.get("AMR_DORFLER", "0.7"))
if AMR_ITERS > 0:
    driven.enable_mesh_refinement(num_iterations=AMR_ITERS, max_DoFs=AMR_MAX_DOFS,
                                  tolerance=AMR_TOL, Dorfler_marking_fraction=AMR_DORFLER,
                                  save_iterations_data=True, save_iterations_mesh=False)
    print(f"Solution-based AMR ENABLED: iters={AMR_ITERS} tol={AMR_TOL:.1e} "
         f"max_DoFs={AMR_MAX_DOFS} Dorfler={AMR_DORFLER} (base mesh unchanged: "
         f"lateral_min={lateral_min_um:.3f}um, i.e. mesh_factor=1.0 always for AMR runs)",
         flush=True)

fs, fe = FSTART_GHZ * 1e9, FEND_GHZ * 1e9
driven.set_freq_values(freq_start=fs, freq_end=fe, freq_step=(fe - fs) / max(N_POINTS - 1, 1))

tp = time.time()
driven.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

# USE_PROM=0 skips patching AdaptiveTol/AdaptiveMaxSamples into the config at
# all, so PALACE falls back to whatever its own Frequency_Driven_Simulation
# default is (per SQuADDS Tutorial 7's setup, likely a plain brute-force sweep
# at every N_POINTS frequency, no ROM) -- added specifically to isolate
# whether AMR's segfault (job 1938861) comes from combining spatial AMR with
# OUR manually-patched PROM config, since that combination has never been
# used elsewhere in this repo (all prior successful AMR runs were on
# capacitance/inductance/eigenmode sims, never a driven+PROM sweep).
USE_PROM = os.environ.get("USE_PROM", "1") != "0"
cfg_path = driven._sim_config
with open(cfg_path) as f:
    cfg = json.load(f)
if USE_PROM:
    cfg["Solver"]["Driven"]["AdaptiveTol"] = ADAPTIVE_TOL
    cfg["Solver"]["Driven"]["AdaptiveMaxSamples"] = ADAPTIVE_MAX_SAMPLES
    with open(cfg_path, "w") as f:
        json.dump(cfg, f, indent=2)
    print(f"PROM: {cfg['Solver']['Driven']}", flush=True)
else:
    print(f"PROM DISABLED (USE_PROM=0): plain sweep, Driven config unchanged: "
         f"{cfg['Solver']['Driven']}", flush=True)

ts = time.time()
driven.run()
solve_s = time.time() - ts
print(f"SOLVE: {solve_s:.1f} s", flush=True)

s = driven.retrieve_data()
if s is None:
    raise SystemExit("retrieve_data() returned None -- the Palace solve failed "
                     "(check the run log for a segfault / SuperLU MemType error).")

# Palace prints its own peak-memory estimate (aggregated correctly across all
# MPI ranks) at the end of out.log -- read that directly rather than trust an
# external `ps --ppid $$` poller in the submit script, which only sees the
# shell's direct children and misses the actual mpirun-spawned PALACE rank
# processes (grandchildren), undercounting real usage for any MPI run.
import re
peak_mem_per_rank_gb = peak_mem_total_gb = None
try:
    log_path = os.path.join(driven._output_data_dir, "out.log")
    with open(log_path) as f:
        log_text = f.read()
    m = re.search(r"Estimated peak per-rank memory usage is:.*?Total ([\d.]+)G", log_text)
    if m:
        peak_mem_per_rank_gb = float(m.group(1))
    m = re.search(r"Estimated peak per-node memory usage is:.*?Total ([\d.]+)G", log_text)
    if m:
        peak_mem_total_gb = float(m.group(1))
    print(f"PALACE-reported peak memory: per-rank-total={peak_mem_per_rank_gb}G "
         f"per-node-total={peak_mem_total_gb}G", flush=True)
except Exception as e:
    print(f"could not parse peak memory from {log_path}: {e}", flush=True)

import numpy as np
freqs = np.asarray(s["freqs"]); S21 = np.asarray(s["S21"])
np.savez(f"{OUTDIR}/results/{TAG}_s21.npz", freq_hz=freqs, S21=S21,
        coupling_space_um=COUPLING_SPACE_UM, coupling_length_um=COUPLING_LENGTH_UM)

from qiskit_metal.analyses.em.transmission_fitting import fit_transmission

# dip_freq_ghz/dip_depth_db are read DIRECTLY off the raw data (argmin of
# |S21|), independent of any Lorentzian fit -- a robustness cross-check.
# Computed BEFORE fitting: a shallow dip (order -1dB, as the first two failed
# pilots showed -- the window missed the true resonance entirely) means
# fit_transmission has nothing but a flat/monotonic trend to fit and its
# curve_fit initial guess can land outside its own parameter bounds
# (ValueError: Initial guess is outside of provided bounds), killing the job
# after the expensive PALACE solve already succeeded. Skip the fit rather than
# let that crash the point -- shallow-dip rows below the threshold get
# Qr/Qc/Qi=None and need a re-run with a shifted/widened window, not a number.
DIP_DEPTH_THRESHOLD_DB = -3.0
s21_db = 20 * np.log10(np.abs(S21))
i_min = int(np.argmin(s21_db))
dip_freq_ghz = float(freqs[i_min] / 1e9)
dip_depth_db = float(s21_db[i_min])

if dip_depth_db > DIP_DEPTH_THRESHOLD_DB:
    print(f"[{TAG}] WARNING: dip_depth_db={dip_depth_db:.3f} > {DIP_DEPTH_THRESHOLD_DB} "
         f"-- window [{FSTART_GHZ},{FEND_GHZ}]GHz likely missed the true resonance "
         f"(dip found only at edge f={dip_freq_ghz}GHz). Skipping fit_transmission.",
         flush=True)
    Qr = Qc = Qi = fr = None
else:
    # fit_transmission returns (dict, plots_list[, raw_params, cov] if full_output) --
    # NOT a bare dict. plot=True returns [fig, ax] (a 3-panel |S21|/phase/IQ figure
    # with the raw data and the fit overlaid) -- save it so a shallow/bad fit (see
    # the module docstring's FREQUENCY WINDOW section) is visible without having
    # to re-fetch and re-plot the raw *_s21.npz by hand.
    fit, plots = fit_transmission(freqs, S21, plot=True)
    Qr, Qc, fr = float(fit["Qr"]), float(fit["Qc"]), float(fit["fr"])
    Qi = 1.0 / (1.0 / Qr - 1.0 / Qc) if Qr != Qc else float("inf")
    if plots:
        plots[0].savefig(f"{OUTDIR}/results/{TAG}_fit.png", dpi=150, bbox_inches="tight")
    fr = fr / 1e9

row = {
    "tag": TAG, "coupling_space_um": COUPLING_SPACE_UM, "coupling_length_um": COUPLING_LENGTH_UM,
    "fr_ghz": fr, "Qr": Qr, "Qc": Qc, "Qi": Qi,
    "dip_freq_ghz": dip_freq_ghz, "dip_depth_db": dip_depth_db,
    "fstart_ghz": FSTART_GHZ, "fend_ghz": FEND_GHZ, "n_points": N_POINTS,
    "adaptive_max_samples": ADAPTIVE_MAX_SAMPLES, "num_cpus": NCPU,
    "solver_tol": SOLVER_TOL, "solver_maxits": SOLVER_MAXITS, "mesh_factor": MESH_FACTOR,
    "amr_iters": AMR_ITERS, "amr_tol": AMR_TOL, "amr_max_dofs": AMR_MAX_DOFS,
    "amr_dorfler": AMR_DORFLER,
    "peak_mem_per_rank_total_gb": peak_mem_per_rank_gb, "peak_mem_per_node_gb": peak_mem_total_gb,
    "solve_s": solve_s, "wall_s": time.time() - t0,
}
with open(f"{OUTDIR}/results/{TAG}.json", "w") as f:
    json.dump(row, f, indent=2)

print(json.dumps(row, indent=2), flush=True)

import glob
for pat in (f"{driven._output_data_dir}/../*.msh", f"{driven._output_data_dir}/paraview"):
    for p in glob.glob(pat):
        os.system(f"rm -rf '{p}'")
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
