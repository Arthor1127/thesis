"""
Electrostatic capacitance-matrix extraction for the N-transmon capacitive-
coupler chain (build_chain.build_chain).

Runs a Palace electrostatic (simPC) capacitance simulation via SQDMetal,
reads the Maxwell matrix straight from Palace's terminal-C.csv (independent of
the SQDMetal post-processing plots that can crash headless), and writes
{OUTDIR}/results/{TAG}.json plus copies of the terminal-*.csv files.

Env:
  N_TRANSMONS      (default 6)
  COUPLER_RES_UM   fine-mesh min size on the couplers (default 3.0)
  TRANSMON_RES_UM  fine-mesh min size on the transmon pads (default 12.0)
  SOLVER_ORDER     (default 2)
  BG_MAX_UM        background max mesh size (default 80)
  MARGIN_X_MM / MARGIN_Y_MM   chip padding beyond the geometry bbox (default 1.2 / 1.2)
  EPS_R            (silicon; SQDMetal 'silicon' -> 11.45)  -- informational only
  PALACE_BIN, NSLOTS, OUTDIR, TAG
No ParaView output (solns_to_save = 0).
"""
import os
import sys
import json
import time
import glob

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

N          = int(os.environ.get("N_TRANSMONS", "6"))
CPL_RES    = float(os.environ.get("COUPLER_RES_UM", "3.0"))
TRM_RES    = float(os.environ.get("TRANSMON_RES_UM", "12.0"))
ORDER      = int(os.environ.get("SOLVER_ORDER", "2"))
BG_MAX_UM  = float(os.environ.get("BG_MAX_UM", "80"))
MARGIN_X   = float(os.environ.get("MARGIN_X_MM", "1.2"))
MARGIN_Y   = float(os.environ.get("MARGIN_Y_MM", "1.2"))
OUTDIR     = os.environ.get("OUTDIR", "./chain_cap")
PALACE_BIN = os.environ["PALACE_BIN"]
NCPU       = int(os.environ.get("NSLOTS", "8"))
TAG        = os.environ.get("TAG", f"chain{N}_c{CPL_RES:g}_o{ORDER}")
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

print(f"[{TAG}] N={N} coupler_res={CPL_RES}um transmon_res={TRM_RES}um "
      f"order={ORDER} bg_max={BG_MAX_UM}um ncpu={NCPU}", flush=True)

t0 = time.time()
from qiskit_metal import designs, draw, Dict            # noqa: F401
from qiskit_metal.qlibrary.core import QComponent        # noqa: F401
# qiskit_metal flips the mpl backend to QtAgg on import; force Agg back BEFORE
# the SQDMetal import or its headless plot calls crash with a QtAgg ImportError.
import matplotlib
matplotlib.use("Agg", force=True)

import build_chain
from build_design import _design_bounds_mm

design = build_chain.build_chain(n_transmons=N)
# size the chip tightly around the geometry (the default margin factor 3 leaves
# ~5 mm of empty ground plane along x for no accuracy gain)
minx, miny, maxx, maxy = _design_bounds_mm(design)
sx, sy = (maxx - minx) + MARGIN_X, (maxy - miny) + MARGIN_Y
cx, cy = 0.5 * (minx + maxx), 0.5 * (miny + maxy)
design.chips.main.size.size_x = f"{sx}mm"
design.chips.main.size.size_y = f"{sy}mm"
design.chips.main.size.center_x = f"{cx}mm"
design.chips.main.size.center_y = f"{cy}mm"
design.rebuild()

comps = list(design.components.keys())
transmons = [c for c in comps if c.startswith("Q")]
couplers  = [c for c in comps if c.startswith("coupler")]
print(f"geometry built ({time.time()-t0:.1f}s)  components={comps}", flush=True)
print(f"chip: {design.chips.main.size.size_x} x {design.chips.main.size.size_y} "
      f"(bbox {maxx-minx:.3f} x {maxy-miny:.3f} mm)", flush=True)

from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation

cap = PALACE_Capacitance_Simulation(
    name=f"cc_{TAG}", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "mesh_refinement": 0,
        "dielectric_material": "silicon",
        "solver_order": ORDER,
        "solver_tol": 1e-8,
        "solver_maxits": 300,
        "solns_to_save": 0,                 # NO ParaView output
        "fillet_resolution": 12,
        "mesh_max": BG_MAX_UM * 1e-6,
        "palace_dir": PALACE_BIN,
        "num_cpus": NCPU,
    })
cap.add_metallic(1)
cap.add_ground_plane()

cap.fine_mesh_components(transmons, min_size=TRM_RES * 1e-6, max_size=BG_MAX_UM * 1e-6,
                         taper_dist_min=10e-6, taper_dist_max=60e-6, metals_only=False)
cap.fine_mesh_components(couplers, min_size=CPL_RES * 1e-6, max_size=BG_MAX_UM * 1e-6,
                         taper_dist_min=5e-6, taper_dist_max=30e-6, metals_only=False)

tp = time.time()
cap.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

ts = time.time()
try:
    cap.run()
except Exception as e:
    print(f"cap.run() post-processing raised (solve output still written): {e}", flush=True)
solve_s = time.time() - ts
print(f"SOLVE: {solve_s:.1f} s", flush=True)

# --- read the Maxwell matrix straight from Palace ---
import pandas as pd
odir = cap._output_data_dir
Craw = pd.read_csv(f"{odir}/terminal-C.csv")
names = [c.strip() for c in Craw.columns[1:]]
C = Craw.to_numpy()[:, 1:].astype(float)          # F, Maxwell 'terminal-C'
try:
    Cm = pd.read_csv(f"{odir}/terminal-Cm.csv").to_numpy()[:, 1:].astype(float)
except Exception:
    Cm = None

# identify ground = the terminal with the largest self-capacitance
gnd = int(np.argmax(np.diag(C)))
n = C.shape[0]
isl = [i for i in range(n) if i != gnd]

# island -> ground and nearest-neighbour couplings (fF)
C_to_gnd_fF   = {names[i]: abs(C[i].sum()) * 1e15 for i in isl}   # row sum = cap to implicit ref
C_selfdiag_fF = {names[i]: C[i, i] * 1e15 for i in isl}
mutual_fF = {}
for a_i in range(n):
    for b_i in range(a_i + 1, n):
        mutual_fF[f"{names[a_i]}--{names[b_i]}"] = abs(C[a_i, b_i]) * 1e15

try:
    sizes = cap.retrieve_simulation_sizes()
except Exception:
    sizes = {}

row = {
    "tag": TAG, "n_transmons": N,
    "coupler_res_um": CPL_RES, "transmon_res_um": TRM_RES,
    "solver_order": ORDER, "bg_max_um": BG_MAX_UM,
    "chip_mm": [sx, sy], "eps_r_sim": 11.45,
    "terminal_names": names, "ground_index": gnd,
    "maxwell_C_matrix_F": C.tolist(),
    "spice_Cm_matrix_F": (Cm.tolist() if Cm is not None else None),
    "C_island_to_ground_fF": C_to_gnd_fF,
    "C_selfdiag_fF": C_selfdiag_fF,
    "C_mutual_fF": mutual_fF,
    "dof": sizes.get("DoF"), "mesh_elements": sizes.get("MeshElements"),
    "solve_s": solve_s, "wall_s": time.time() - t0,
}
with open(f"{OUTDIR}/results/{TAG}.json", "w") as f:
    json.dump(row, f, indent=2)
for src in ("terminal-C.csv", "terminal-Cm.csv", "terminal-Cinv.csv"):
    p = f"{odir}/{src}"
    if os.path.exists(p):
        os.system(f"cp '{p}' '{OUTDIR}/results/{TAG}_{src}'")

print(json.dumps({k: v for k, v in row.items()
                  if k not in ("maxwell_C_matrix_F", "spice_Cm_matrix_F")}, indent=2), flush=True)

# drop the heavy mesh; keep CSVs + json
for pat in (f"{odir}/../*.msh", f"{odir}/paraview"):
    for p in glob.glob(pat):
        os.system(f"rm -rf '{p}'")
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
