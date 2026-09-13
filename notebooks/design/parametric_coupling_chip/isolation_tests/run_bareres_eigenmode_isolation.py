"""Isolation test for the eigenmode Qc-from-lossy-port method (see
run_qc_eigenmode.py's docstring): a BARE quarter-wave CPW resonator --
ShortToGround at one end, a SINGLE lumped port of resistance R at the other,
NO CoupledLineTee, NO launchpads, NO snake routing -- with a clean, standard
ANALYTIC prediction to check against.

WHY: the first real-design eigenmode run (cs14_cl300_eig) found Palace's own
eig.csv "loaded_Q" (=228) and its own port-Q.csv-derived Qc (=fr/kappa_ext,
summed over both ports =13562) disagree with EACH OTHER by ~60x for the same
mode, even though the port-Q.csv number happened to land within ~13.5% of
the driven-S21 fit (15684). That is not enough to trust the method -- need a
case with a KNOWN right answer and zero coupler-geometry ambiguity.

ANALYTIC PREDICTION: for a shorted quarter-wave transmission-line resonator
(impedance Z0) loaded by a resistor R at its open end, standard cavity
perturbation theory (Pozar, Microwave Engineering) gives the effective
capacitance at the open-end voltage antinode as C_eff = pi/(4*omega_0*Z0),
so a shunt conductance G=1/R there gives:
    Qc = omega_0 * C_eff / G = (pi/4) * (R/Z0)
This is the SAME formula family already used elsewhere in this repo's own
design_helper.py (compute_resonator_Qc_from_Cc), just for a resistive
(not capacitive) load -- no coupler geometry, no coupling capacitance to
estimate, nothing to get wrong on the analytic side.

R is chosen per-run to land Qc_analytic near our real design's regime
(~1e4-1e5), since that's exactly where the driven-sweep noise-floor problem
showed up -- a validation at, say, Qc~10 would prove nothing about whether
the method works at Qc~15000.

Env: PORT_R_OHM, PALACE_BIN, [NSLOTS=16], [OUTDIR=./out_bareres_eig], [TAG],
     [STARTING_FREQ_GHZ=7.5], [NUM_FREQS=2], [SOLVER_ORDER=1], [SOLVER_TOL=1e-6],
     [SOLVER_MAXITS=100]
Writes: {OUTDIR}/results/{TAG}.json
Cluster job only -- run via submit_bareres_eigenmode.sh.
"""
import os
import sys
import json
import time
import math

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PORT_R_OHM   = float(os.environ["PORT_R_OHM"])
OUTDIR       = os.environ.get("OUTDIR", "./out_bareres_eig")
PALACE_BIN   = os.environ["PALACE_BIN"]
NCPU         = int(os.environ.get("NSLOTS", "16"))
TAG          = os.environ.get("TAG", f"R{PORT_R_OHM:g}")
STARTING_FREQ_GHZ = float(os.environ.get("STARTING_FREQ_GHZ", "7.5"))
NUM_FREQS    = int(os.environ.get("NUM_FREQS", "2"))
SOLVER_ORDER = int(os.environ.get("SOLVER_ORDER", "1"))
SOLVER_TOL   = float(os.environ.get("SOLVER_TOL", "1e-6"))
SOLVER_MAXITS = int(os.environ.get("SOLVER_MAXITS", "100"))
# SQDMetal hardcodes Solver.Linear.Type="SuperLU" / KSPType="FGMRES" for
# eigenmode. PalaceForCQED (the trusted reference) uses "Default"/"GMRES" --
# Palace's "Default" picks a proper AMS/multigrid preconditioner for the
# curl-curl operator, which SuperLU-as-configured here does NOT provide (job
# 1938974's order-4 runs showed bare unpreconditioned FGMRES crawling for
# ~8 h). Set LINEAR_TYPE / KSP_TYPE to patch the generated config.json after
# prepare_simulation(); leave empty to keep SQDMetal's default.
LINEAR_TYPE  = os.environ.get("LINEAR_TYPE", "").strip()
KSP_TYPE     = os.environ.get("KSP_TYPE", "").strip()
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

CPW_WIDTH_UM = 20.0
CPW_GAP_UM = 11.5
Z0_OHM = 49.961   # DesignHelper.compute_cpw_geometric_parameters at this cross-section
                  # (readout_resonators.py's own validated reference value)
TARGET_RES_LENGTH_MM = 3.77259   # same lambda/4-at-8GHz length as the real design
                                 # (readout_resonators.py), for a directly comparable regime

Qc_analytic = (math.pi / 4.0) * (PORT_R_OHM / Z0_OHM)
print(f"[{TAG}] PORT_R_OHM={PORT_R_OHM:g} Z0={Z0_OHM} -> Qc_analytic={Qc_analytic:.2f} "
     f"ncpu={NCPU} solver_order={SOLVER_ORDER} solver_tol={SOLVER_TOL:g} "
     f"solver_maxits={SOLVER_MAXITS}", flush=True)

t0 = time.time()
from qiskit_metal import designs, Dict
import matplotlib
matplotlib.use("Agg", force=True)

from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder

half_span_um = TARGET_RES_LENGTH_MM * 1e3 / 2.0
chip_size_x_mm = (TARGET_RES_LENGTH_MM * 1e3 + 800.0) / 1e3   # +400um margin each side

design = designs.DesignPlanar({}, overwrite_enabled=True)
design.chips.main.size.size_x = f'{chip_size_x_mm}mm'
design.chips.main.size.size_y = '2mm'
design.chips.main.size.size_z = '500um'
design.chips.main.size.center_x = '0mm'
design.chips.main.size.center_y = '0mm'
design.variables['cpw_width'] = f'{CPW_WIDTH_UM}um'
design.variables['cpw_gap'] = f'{CPW_GAP_UM}um'

# OTG1 (open end, where the lumped R port attaches) -- ShortToGround (the
# real, physical short) at the other end. Same OpenToGround-as-a-
# RoutePathfinder-anchor trick as check_line_impedance.py's line_only mode:
# the port attaches to the ROUTE's own pin, not to OTG1 itself, so OTG1's
# only job is giving RoutePathfinder a valid pin to route from.
x1 = f'-{half_span_um}um'
x2 = f'{half_span_um}um'
OpenToGround(design, 'OTG1', options=dict(
    chip='main', pos_x=x1, pos_y='0um', orientation='180', width='cpw_width'))
ShortToGround(design, 'STG', options=dict(
    chip='main', pos_x=x2, pos_y='0um', orientation='0', width='cpw_width'))
RoutePathfinder(design, 'RES', options=dict(
    chip='main', trace_width='cpw_width', trace_gap='cpw_gap', fillet='50um',
    pin_inputs=Dict(start_pin=Dict(component='OTG1', pin='open'),
                    end_pin=Dict(component='STG', pin='short'))))
design.rebuild()
print(f"geometry built ({time.time()-t0:.1f}s): OTG1--RES--STG, length~"
     f"{TARGET_RES_LENGTH_MM}mm, chip={design.chips.main.size.size_x} x "
     f"{design.chips.main.size.size_y}", flush=True)

from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation

eigen = PALACE_Eigenmode_Simulation(
    name=f"bareres_{TAG}", metal_design=design, sim_parent_directory=OUTDIR,
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
# ONLY ONE port -- the short at the other end is a real PEC short (lossless
# by construction), not a second lossy channel to combine.
eigen.create_port_CPW_on_Route("RES", pin_name="start", len_launch=30e-6,
                               impedance_R=PORT_R_OHM)

lateral_min_um = 0.5 * min(CPW_WIDTH_UM, CPW_GAP_UM)
eigen.fine_mesh_components(["RES", "OTG1", "STG"], min_size=lateral_min_um * 1e-6,
                           max_size=150e-6, taper_dist_min=16e-6, taper_dist_max=40e-6,
                           metals_only=False)
print(f"mesh: lateral_min={lateral_min_um:.3f}um", flush=True)

tp = time.time()
eigen.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

if LINEAR_TYPE or KSP_TYPE:
    cfg_path = eigen._sim_config
    with open(cfg_path) as f:
        cfg = json.load(f)
    if LINEAR_TYPE:
        cfg["Solver"]["Linear"]["Type"] = LINEAR_TYPE
    if KSP_TYPE:
        cfg["Solver"]["Linear"]["KSPType"] = KSP_TYPE
    with open(cfg_path, "w") as f:
        json.dump(cfg, f, indent=2)
    print(f"PATCHED Solver.Linear -> {cfg['Solver']['Linear']}", flush=True)

ts = time.time()
eigen.run()
solve_s = time.time() - ts
print(f"SOLVE: {solve_s:.1f} s", flush=True)

eig_csv_path = eigen._output_data_dir + "/eig.csv"
if not os.path.exists(eig_csv_path):
    raise SystemExit(f"FAILED: {eig_csv_path} was never written (job likely killed mid-run).")

port_epr = eigen.retrieve_mode_port_EPR()
freqs_hz = port_epr["eigenfrequencies"]
loaded_Q = port_epr["loaded_Q"]
kappa_hz = port_epr["kappa"] * 1e9   # port-Q.csv is in GHz despite the field name -- see
                                     # run_qc_eigenmode.py's fix/comment for the same bug.

import numpy as np
n_modes = len(freqs_hz)
if n_modes == 0:
    raise SystemExit(f"FAILED: zero converged modes (starting_freq={STARTING_FREQ_GHZ}GHz).")

modes = []
for i in range(n_modes):
    k = float(kappa_hz[i]) if np.isfinite(kappa_hz[i]) else None
    f_ghz = float(freqs_hz[i]) / 1e9
    Qc_i = f_ghz * 1e9 / k if (k is not None and k > 0) else None
    ratio = (loaded_Q[i] / Qc_i) if (Qc_i and loaded_Q[i] > 0) else None
    modes.append({"m": i + 1, "f_ghz": f_ghz, "loaded_Q": float(loaded_Q[i]),
                 "kappa_hz": k, "Qc_from_kappa": Qc_i, "loaded_Q_over_Qc_kappa": ratio})
    print(f"  mode {i+1}: f={f_ghz:.5f}GHz loaded_Q={loaded_Q[i]:.4e} kappa={k}Hz "
         f"Qc_from_kappa={Qc_i} loaded_Q/Qc_kappa={ratio}", flush=True)

# Primary mode = the one nearest starting_freq (only one port -> no
# largest-kappa ambiguity needed, but multiple modes can still converge).
primary = min(modes, key=lambda m: abs(m["f_ghz"] - STARTING_FREQ_GHZ))

row = {
    "tag": TAG, "port_R_ohm": PORT_R_OHM, "Z0_ohm": Z0_OHM, "Qc_analytic": Qc_analytic,
    "fr_ghz": primary["f_ghz"], "loaded_Q": primary["loaded_Q"],
    "kappa_hz": primary["kappa_hz"], "Qc_from_kappa": primary["Qc_from_kappa"],
    "loaded_Q_over_Qc_kappa": primary["loaded_Q_over_Qc_kappa"],
    "Qc_from_kappa_over_analytic": (primary["Qc_from_kappa"] / Qc_analytic
                                    if primary["Qc_from_kappa"] else None),
    "loaded_Q_over_analytic": primary["loaded_Q"] / Qc_analytic,
    "modes": modes,
    "starting_freq_ghz": STARTING_FREQ_GHZ, "number_of_freqs": NUM_FREQS, "num_cpus": NCPU,
    "solver_order": SOLVER_ORDER, "solver_tol": SOLVER_TOL, "solver_maxits": SOLVER_MAXITS,
    "linear_type": LINEAR_TYPE or "SuperLU", "ksp_type": KSP_TYPE or "FGMRES",
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
