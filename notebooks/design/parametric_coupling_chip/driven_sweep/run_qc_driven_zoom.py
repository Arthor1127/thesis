"""Phase D of the Qc investigation: PalaceForCQED-style *narrow-window* driven
S21 with a complex hanger-notch fit, for the design points where the eigenmode
method saturates (Qc above its ~85k-1e-6-tol / higher-order ceiling -- see
../QC_STATUS.md, ../HANDOFF.md).

Difference from run_qc_sim.py (the wide-window production driver):
  * The frequency window is a THIN band centred on a known f_r -- pass
    FR_SEED_GHZ (from the eigenmode result's `fr_ghz`, or a fast order-1
    eigenmode pre-pass). Default half-window 10 MHz. This is the recipe from
    PalaceForCQED's examples/run.jl narrow scans (7-10 MHz windows, 20 kHz
    step) -- the previous attempts here used a 750 MHz window that defeated
    PROM's adaptive sampler (see run_qc_sim.py's FREQUENCY WINDOW note and the
    palace-driven-sweep-resonance-finding memory).
  * The fit is fit_qc_notch.fit_notch (complex S21, diameter-corrected hanger
    model, cable-delay + environment pre-removal), with qiskit_metal's
    fit_transmission also run for comparison.

Env: COUPLING_SPACE_UM, COUPLING_LENGTH_UM, PALACE_BIN, FR_SEED_GHZ,
     [NSLOTS=16], [OUTDIR=./out_qc_zoom], [TAG],
     [HALFWIN_MHZ=10], [NPTS=2001], [ADAPTIVE_TOL=1e-3],
     [ADAPTIVE_MAX_SAMPLES=400], [SOLVER_ORDER=1], [SOLVER_TOL=1e-8],
     [SOLVER_MAXITS=300], [MESH_FACTOR=0.75]
Writes: {OUTDIR}/results/{TAG}.json, {TAG}_s21.npz, {TAG}_fit.png

Cluster job only -- run via submit_qc_driven_zoom.sh.
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root

COUPLING_SPACE_UM  = float(os.environ["COUPLING_SPACE_UM"])
COUPLING_LENGTH_UM = float(os.environ["COUPLING_LENGTH_UM"])
OUTDIR       = os.environ.get("OUTDIR", "./out_qc_zoom")
PALACE_BIN   = os.environ["PALACE_BIN"]
NCPU         = int(os.environ.get("NSLOTS", "16"))
TAG          = os.environ.get("TAG", f"cs{COUPLING_SPACE_UM:g}_cl{COUPLING_LENGTH_UM:g}")
FR_SEED_GHZ  = float(os.environ["FR_SEED_GHZ"])
HALFWIN_MHZ  = float(os.environ.get("HALFWIN_MHZ", "10"))
N_POINTS     = int(os.environ.get("NPTS", "2001"))
ADAPTIVE_TOL = float(os.environ.get("ADAPTIVE_TOL", "1e-3"))
ADAPTIVE_MAX_SAMPLES = int(os.environ.get("ADAPTIVE_MAX_SAMPLES", "400"))
SOLVER_ORDER = int(os.environ.get("SOLVER_ORDER", "1"))
SOLVER_TOL   = float(os.environ.get("SOLVER_TOL", "1e-8"))
SOLVER_MAXITS = int(os.environ.get("SOLVER_MAXITS", "300"))
MESH_FACTOR  = float(os.environ.get("MESH_FACTOR", "0.75"))
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

FSTART_GHZ = FR_SEED_GHZ - HALFWIN_MHZ * 1e-3
FEND_GHZ   = FR_SEED_GHZ + HALFWIN_MHZ * 1e-3
step_khz = (FEND_GHZ - FSTART_GHZ) * 1e6 / max(N_POINTS - 1, 1)
print(f"[{TAG}] cs={COUPLING_SPACE_UM}um cl={COUPLING_LENGTH_UM}um  fr_seed={FR_SEED_GHZ}GHz "
      f"window [{FSTART_GHZ:.6f},{FEND_GHZ:.6f}]GHz ({2*HALFWIN_MHZ:g}MHz, {step_khz:.2f}kHz step, "
      f"{N_POINTS}pts)  order={SOLVER_ORDER} tol={SOLVER_TOL:g}  ncpu={NCPU}", flush=True)

t0 = time.time()
from qiskit_metal import designs  # noqa: F401  (forces the QtAgg->Agg reorder below)
import matplotlib
matplotlib.use("Agg", force=True)

from readout_resonators import build_design

design = build_design(COUPLING_SPACE_UM, COUPLING_LENGTH_UM)
print(f"geometry built ({time.time()-t0:.1f}s)", flush=True)

from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation

driven = PALACE_Driven_Simulation(
    name=f"qcz_{TAG}", metal_design=design, sim_parent_directory=OUTDIR,
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

# Same mesh recipe as run_qc_sim.py: the coupling gap itself sets the fine-mesh
# floor (palace-mesh-resolution-and-amr).
lateral_min_um = 0.5 * MESH_FACTOR * min(20.0, 11.5, COUPLING_SPACE_UM)
res_min_um = max(lateral_min_um, 1.0)
bg_max_um = 200.0
driven.fine_mesh_components(["TL1", "TL2", "LP1", "LP2"], min_size=8.0e-6, max_size=bg_max_um*1e-6,
                            taper_dist_min=16e-6, taper_dist_max=40e-6, metals_only=False)
driven.fine_mesh_components(["readout_1", "CLT", "res_short"], min_size=res_min_um*1e-6,
                            max_size=bg_max_um*1e-6, taper_dist_min=30e-6, taper_dist_max=100e-6,
                            metals_only=False)
print(f"mesh: lateral_min={lateral_min_um:.3f}um res_min={res_min_um:.3f}um", flush=True)

fs, fe = FSTART_GHZ * 1e9, FEND_GHZ * 1e9
driven.set_freq_values(freq_start=fs, freq_end=fe, freq_step=(fe - fs) / max(N_POINTS - 1, 1))

tp = time.time()
driven.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

cfg_path = driven._sim_config
with open(cfg_path) as f:
    cfg = json.load(f)
cfg["Solver"]["Driven"]["AdaptiveTol"] = ADAPTIVE_TOL
cfg["Solver"]["Driven"]["AdaptiveMaxSamples"] = ADAPTIVE_MAX_SAMPLES
with open(cfg_path, "w") as f:
    json.dump(cfg, f, indent=2)
print(f"PROM: {cfg['Solver']['Driven']}", flush=True)

ts = time.time()
driven.run()
solve_s = time.time() - ts
print(f"SOLVE: {solve_s:.1f} s", flush=True)

s = driven.retrieve_data()
if s is None:
    raise SystemExit("retrieve_data() returned None -- the Palace solve failed.")

import re
peak_mem_per_rank_gb = peak_mem_total_gb = None
try:
    with open(os.path.join(driven._output_data_dir, "out.log")) as f:
        log_text = f.read()
    m = re.search(r"Estimated peak per-rank memory usage is:.*?Total ([\d.]+)G", log_text)
    if m:
        peak_mem_per_rank_gb = float(m.group(1))
    m = re.search(r"Estimated peak per-node memory usage is:.*?Total ([\d.]+)G", log_text)
    if m:
        peak_mem_total_gb = float(m.group(1))
except Exception as e:
    print(f"peak-mem parse failed: {e}", flush=True)

import numpy as np
freqs = np.asarray(s["freqs"]); S21 = np.asarray(s["S21"])
np.savez(f"{OUTDIR}/results/{TAG}_s21.npz", freq_hz=freqs, S21=S21,
         coupling_space_um=COUPLING_SPACE_UM, coupling_length_um=COUPLING_LENGTH_UM,
         fr_seed_ghz=FR_SEED_GHZ)

s21_db = 20 * np.log10(np.abs(S21))
i_min = int(np.argmin(s21_db))
dip_freq_ghz = float(freqs[i_min] / 1e9)
dip_depth_db = float(s21_db[i_min] - np.median(s21_db))   # relative to baseline
edge_miss = bool(i_min < 3 or i_min > len(freqs) - 4)

# --- complex notch fit ---
from fit_qc_notch import fit_notch
nf = None
try:
    nf = fit_notch(freqs, S21)
    print(f"notch fit: fr={nf['fr_hz']/1e9:.6f}GHz Ql={nf['Ql']:.0f} Qc={nf['Qc']:.0f} "
          f"Qi={nf['Qi']:.0f} phi={nf['phi']:+.3f} depth={nf['depth_db']:+.2f}dB "
          f"rmse={nf['rmse']:.2e} {'OK' if nf['success'] else 'FAIL'}", flush=True)
except Exception as e:
    print(f"notch fit raised: {e}", flush=True)

# --- qiskit_metal fit for comparison ---
qm_fr = qm_Qr = qm_Qc = qm_Qi = None
try:
    from qiskit_metal.analyses.em.transmission_fitting import fit_transmission
    if not edge_miss and dip_depth_db < -0.3:
        fit, plots = fit_transmission(freqs, S21, plot=False)
        qm_Qr, qm_Qc, qm_fr = float(fit["Qr"]), float(fit["Qc"]), float(fit["fr"]) / 1e9
        qm_Qi = 1.0 / (1.0 / qm_Qr - 1.0 / qm_Qc) if qm_Qr != qm_Qc else float("inf")
except Exception as e:
    print(f"qiskit_metal fit_transmission raised: {e}", flush=True)

# --- diagnostic plot ---
try:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    fghz = freqs / 1e9
    ax[0].plot(fghz, np.abs(S21), ".", ms=3, label="PALACE")
    ax[1].plot(fghz, np.unwrap(np.angle(S21)), ".", ms=3)
    ax[2].plot(S21.real, S21.imag, ".", ms=3)
    if nf is not None:
        ax[0].plot(fghz, np.abs(nf["model"]), "-", lw=1, label="notch fit")
        ax[1].plot(fghz, np.unwrap(np.angle(nf["model"])), "-", lw=1)
        ax[2].plot(nf["model"].real, nf["model"].imag, "-", lw=1)
    ax[0].set(xlabel="f (GHz)", ylabel="|S21|", title=f"{TAG}  Qc(notch)="
              f"{nf['Qc']:.3g}" if nf else TAG); ax[0].legend()
    ax[1].set(xlabel="f (GHz)", ylabel="arg S21 (rad)")
    ax[2].set(xlabel="Re S21", ylabel="Im S21"); ax[2].set_aspect("equal", "box")
    fig.tight_layout()
    fig.savefig(f"{OUTDIR}/results/{TAG}_fit.png", dpi=140, bbox_inches="tight")
except Exception as e:
    print(f"plot failed: {e}", flush=True)

row = {
    "tag": TAG, "coupling_space_um": COUPLING_SPACE_UM, "coupling_length_um": COUPLING_LENGTH_UM,
    "fr_seed_ghz": FR_SEED_GHZ, "halfwin_mhz": HALFWIN_MHZ, "n_points": N_POINTS,
    "step_khz": step_khz,
    "dip_freq_ghz": dip_freq_ghz, "dip_depth_db": dip_depth_db, "edge_miss": edge_miss,
    "notch_fr_ghz": (nf["fr_hz"] / 1e9) if nf else None,
    "notch_Ql": nf["Ql"] if nf else None,
    "notch_Qc": nf["Qc"] if nf else None,
    "notch_Qi": nf["Qi"] if nf else None,
    "notch_phi": nf["phi"] if nf else None,
    "notch_depth_db": nf["depth_db"] if nf else None,
    "notch_rmse": nf["rmse"] if nf else None,
    "notch_success": nf["success"] if nf else None,
    "qm_fr_ghz": qm_fr, "qm_Qr": qm_Qr, "qm_Qc": qm_Qc, "qm_Qi": qm_Qi,
    "adaptive_max_samples": ADAPTIVE_MAX_SAMPLES, "num_cpus": NCPU,
    "solver_order": SOLVER_ORDER, "solver_tol": SOLVER_TOL, "solver_maxits": SOLVER_MAXITS,
    "mesh_factor": MESH_FACTOR,
    "peak_mem_per_rank_total_gb": peak_mem_per_rank_gb, "peak_mem_per_node_gb": peak_mem_total_gb,
    "solve_s": solve_s, "wall_s": time.time() - t0,
}
with open(f"{OUTDIR}/results/{TAG}.json", "w") as f:
    json.dump(row, f, indent=2)
print(json.dumps({k: v for k, v in row.items() if k != "model"}, indent=2), flush=True)

import glob
for pat in (f"{driven._output_data_dir}/../*.msh", f"{driven._output_data_dir}/paraview"):
    for p in glob.glob(pat):
        os.system(f"rm -rf '{p}'")
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
