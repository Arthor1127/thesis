"""Eigenmode Qc extraction for the two-port hanger isolation test (see
hanger_isolation_design.py's docstring for the full motivation): a
side-coupled (not directly loaded) quarter-wave resonator, coupled to a real
2-port 50-Ohm through-line via an idealized lumped capacitor of exactly known
value, with the standard hanger analytic formula
Qc = 2*C_eff/(Z0*omega_r*Cc^2) as ground truth.

This exercises the actual TWO-PORT port-Q.csv kappa-summing path that the
1-port bare-resonator isolation test (run_bareres_eigenmode_isolation.py)
structurally cannot -- that test validated the eigenmode method's numerical
precision vs Qc; this one checks whether the real design's ~60x
loaded_Q-vs-kappa disagreement is a "two lossy ports" issue in general, or
specific to the real design's extra complexity (nearby competing mode
~250MHz away, coupler geometry, snake routing, etc. -- none of which this
design has).

Env: QC_TARGET (300|3000|30000, or any float -- solved via Cc_for_Qc),
     PALACE_BIN, [NSLOTS=16], [OUTDIR=./out_hanger_eig], [TAG],
     [STARTING_FREQ_GHZ=7.3], [NUM_FREQS=2], [SOLVER_TOL=1e-6], [SOLVER_MAXITS=100]
Writes: {OUTDIR}/results/{TAG}.json
Cluster job only.
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

QC_TARGET    = float(os.environ["QC_TARGET"])
OUTDIR       = os.environ.get("OUTDIR", "./out_hanger_eig")
PALACE_BIN   = os.environ["PALACE_BIN"]
NCPU         = int(os.environ.get("NSLOTS", "16"))
TAG          = os.environ.get("TAG", f"hanger_eig_Qc{QC_TARGET:g}")
STARTING_FREQ_GHZ = float(os.environ.get("STARTING_FREQ_GHZ", "7.3"))
NUM_FREQS    = int(os.environ.get("NUM_FREQS", "2"))
SOLVER_TOL   = float(os.environ.get("SOLVER_TOL", "1e-6"))
SOLVER_MAXITS = int(os.environ.get("SOLVER_MAXITS", "100"))
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

t0 = time.time()
from qiskit_metal import designs  # noqa: F401
import matplotlib
matplotlib.use("Agg", force=True)

import hanger_isolation_design as hid

Cc_F = hid.Cc_for_Qc(QC_TARGET)
print(f"[{TAG}] QC_TARGET={QC_TARGET:g} -> Cc={Cc_F*1e15:.4f}fF  ncpu={NCPU}", flush=True)

design = hid.build_design(Cc_F)
print(f"geometry built ({time.time()-t0:.1f}s)  chip={design.chips.main.size.size_x} "
      f"x {design.chips.main.size.size_y}", flush=True)

from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation

eigen = PALACE_Eigenmode_Simulation(
    name=f"heig_{TAG}", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "dielectric_material": "silicon",
        "starting_freq": STARTING_FREQ_GHZ * 1e9, "number_of_freqs": NUM_FREQS,
        "solns_to_save": NUM_FREQS,
        "solver_order": 1, "solver_tol": SOLVER_TOL, "solver_maxits": SOLVER_MAXITS,
        "palace_dir": PALACE_BIN, "num_cpus": NCPU,
    })
eigen.add_metallic(1)
eigen.add_ground_plane()
eigen.create_port_CPW_on_Route("TL1", pin_name="start", len_launch=30e-6, impedance_R=50)
eigen.create_port_CPW_on_Route("TL2", pin_name="end", len_launch=30e-6, impedance_R=50)
eigen.create_port_2_conds("JT", "second_end", "RES", "start",
                          rect_width=20e-6, impedance_R=0, impedance_C=Cc_F)

lateral_min_um = 0.5 * min(hid.CPW_WIDTH_UM, hid.CPW_GAP_UM, hid.COUPLING_GAP_UM)
eigen.fine_mesh_components(["TL1", "TL2", "JT", "OTG1", "OTG2"], min_size=lateral_min_um * 1e-6,
                           max_size=150e-6, taper_dist_min=16e-6, taper_dist_max=40e-6,
                           metals_only=False)
eigen.fine_mesh_components(["RES", "OTGr", "STG"], min_size=lateral_min_um * 1e-6,
                           max_size=150e-6, taper_dist_min=30e-6, taper_dist_max=100e-6,
                           metals_only=False)
print(f"mesh: lateral_min={lateral_min_um:.3f}um", flush=True)

tp = time.time()
eigen.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

ts = time.time()
eigen.run()
solve_s = time.time() - ts
print(f"SOLVE: {solve_s:.1f} s", flush=True)

eig_csv_path = eigen._output_data_dir + "/eig.csv"
if not os.path.exists(eig_csv_path):
    raise SystemExit(f"FAILED: {eig_csv_path} was never written.")

port_epr = eigen.retrieve_mode_port_EPR()
freqs_hz = port_epr["eigenfrequencies"]
loaded_Q = port_epr["loaded_Q"]
kappa_hz = port_epr["kappa"] * 1e9   # port-Q.csv is in GHz despite the field name

import numpy as np
n_modes = len(freqs_hz)
if n_modes == 0:
    raise SystemExit(f"FAILED: zero converged modes.")

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

primary = min(modes, key=lambda m: abs(m["f_ghz"] - STARTING_FREQ_GHZ))

row = {
    "tag": TAG, "qc_target": QC_TARGET, "Cc_fF": Cc_F * 1e15, "Z0_ohm": hid.Z0_OHM,
    "fr_ghz": primary["f_ghz"], "loaded_Q": primary["loaded_Q"],
    "kappa_hz": primary["kappa_hz"], "Qc_from_kappa": primary["Qc_from_kappa"],
    "loaded_Q_over_Qc_kappa": primary["loaded_Q_over_Qc_kappa"],
    "Qc_from_kappa_over_analytic": (primary["Qc_from_kappa"] / QC_TARGET
                                    if primary["Qc_from_kappa"] else None),
    "modes": modes,
    "starting_freq_ghz": STARTING_FREQ_GHZ, "number_of_freqs": NUM_FREQS, "num_cpus": NCPU,
    "solver_tol": SOLVER_TOL, "solver_maxits": SOLVER_MAXITS,
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
