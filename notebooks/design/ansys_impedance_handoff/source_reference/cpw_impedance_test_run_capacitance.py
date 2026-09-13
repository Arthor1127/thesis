"""Electrostatic capacitance-per-unit-length of a straight CPW segment
(trace <-> ground Maxwell matrix), one (cpw_width, cpw_gap[, thickness]) point.

Env: CPW_WIDTH_UM, CPW_GAP_UM, [LENGTH_UM=3000], PALACE_BIN, [NSLOTS=8],
     [OUTDIR=./out_cap], [TAG], [SOLVER_ORDER=1], [BG_MAX_UM=80],
     [FULL3D=0], [METAL_THICKNESS_NM=100] (only used if FULL3D=1),
     [AMR=1], [AMR_ITERS=4], [AMR_TOL=1e-2], [AMR_MAX_DOFS=3000000], [AMR_DORFLER=0.6]
Writes: {OUTDIR}/results/{TAG}.json

Cluster job only -- run via submit_capacitance.sh. Do not run this locally
(see ../../../directives -- PALACE/FEM solves always go to the CNEA cluster).

Mesh sizing: the base fine-mesh element size is HALF the smallest lateral
feature in the design (min(cpw_width, cpw_gap)) -- e.g. 5um for a 10um gap.
For FULL3D (a real metal thickness, not the default zero-thickness sheet),
GMSH needs a real element across that thickness too; per
scripts/python/island_cap/run_cap_island.py's finding, forcing res_min down to
"half the thickness" for a 20-300nm film would demand few-to-sub-nm elements
over a multi-mm trace (unsolvable), so the through-thickness floor instead
follows that script's convention: res_min = max(half-smallest-lateral-feature,
metal_thickness, 0.15um) -- effectively still the lateral rule in every case
swept here, since these thicknesses are all << the lateral mesh size anyway.
AMR (enable_mesh_refinement) then refines from that starting point toward a
target relative error, rather than trusting the fixed base mesh alone.
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CPW_WIDTH_UM = float(os.environ["CPW_WIDTH_UM"])
CPW_GAP_UM   = float(os.environ["CPW_GAP_UM"])
LENGTH_UM    = float(os.environ.get("LENGTH_UM", "3000"))
OUTDIR       = os.environ.get("OUTDIR", "./out_cap")
PALACE_BIN   = os.environ["PALACE_BIN"]
NCPU         = int(os.environ.get("NSLOTS", "8"))
TAG          = os.environ.get("TAG", f"w{CPW_WIDTH_UM:g}_g{CPW_GAP_UM:g}")
BG_MAX_UM    = float(os.environ.get("BG_MAX_UM", "80"))
FULL3D       = os.environ.get("FULL3D", "0") == "1"
METAL_THICKNESS_NM = float(os.environ.get("METAL_THICKNESS_NM", "100"))
AMR          = os.environ.get("AMR", "1") == "1"
AMR_ITERS    = int(os.environ.get("AMR_ITERS", "4"))
AMR_TOL      = float(os.environ.get("AMR_TOL", "1e-2"))
AMR_MAX_DOFS = int(os.environ.get("AMR_MAX_DOFS", "3000000"))
AMR_DORFLER  = float(os.environ.get("AMR_DORFLER", "0.6"))
DIELECTRIC_MATERIAL = os.environ.get("DIELECTRIC_MATERIAL", "silicon")
if DIELECTRIC_MATERIAL.lower() == "vacuum":
    DIELECTRIC_MATERIAL = ""  # Material("") -> er=1, per=1, tan_delta=0 (SQDMetal's own default branch)
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

print(f"[{TAG}] cpw_width={CPW_WIDTH_UM}um cpw_gap={CPW_GAP_UM}um "
      f"length={LENGTH_UM}um ncpu={NCPU} FULL3D={FULL3D} "
      f"thickness={METAL_THICKNESS_NM if FULL3D else 'n/a (zero-thickness)'}nm "
      f"AMR={AMR}", flush=True)

t0 = time.time()
from qiskit_metal import designs  # noqa: F401  (forces the QtAgg->Agg reorder below)
import matplotlib
matplotlib.use("Agg", force=True)

from build_design import build_design, OTG_NAME, STG_NAME, ROUTE_NAME

# far_end='open': the trace MUST be a floating conductor for a capacitance
# extraction. far_end='short' (run_inductance.py's geometry) galvanically
# shorts trace to ground, which would collapse this into a single-conductor,
# undefined-capacitance problem instead -- see build_design.build_design's
# docstring.
design = build_design(CPW_WIDTH_UM, CPW_GAP_UM, length_um=LENGTH_UM, far_end='open')
print(f"geometry built ({time.time()-t0:.1f}s)  chip={design.chips.main.size.size_x} "
      f"x {design.chips.main.size.size_y}", flush=True)

from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation

cap = PALACE_Capacitance_Simulation(
    name=f"cc_{TAG}", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "dielectric_material": DIELECTRIC_MATERIAL,
        "solver_order": int(os.environ.get("SOLVER_ORDER", "1")),
        "solver_tol": 1e-8,
        "solver_maxits": 200,
        "solns_to_save": 0,
        "fillet_resolution": 10,
        "mesh_max": BG_MAX_UM * 1e-6,
        "palace_dir": PALACE_BIN,
        "num_cpus": NCPU,
    })
cap.add_metallic(1)
cap.add_ground_plane()

if FULL3D:
    cap.enforce_full_3D_simulation(metal_thickness=METAL_THICKNESS_NM * 1e-9,
                                   substrate_trenching=0.0)

# base fine mesh: half the smallest lateral feature (see module docstring)
lateral_min_um = 0.5 * min(CPW_WIDTH_UM, CPW_GAP_UM)
res_min_um = max(lateral_min_um, METAL_THICKNESS_NM * 1e-3 if FULL3D else 0.0, 0.15)
res_min = res_min_um * 1e-6
taper_max = max(4.0 * CPW_GAP_UM, 20.0) * 1e-6
# The open/short termination gaps must be in this fine-mesh list too
# (sqdmetal-termination-gap-meshing lesson) or GMSH under-resolves the notch
# and an "open" end shorts to ground.
cap.fine_mesh_components([ROUTE_NAME, OTG_NAME, STG_NAME], min_size=res_min,
                         max_size=BG_MAX_UM * 1e-6,
                         taper_dist_min=1.5 * CPW_GAP_UM * 1e-6, taper_dist_max=taper_max,
                         metals_only=False)
print(f"mesh: lateral_min={lateral_min_um:.3f}um  res_min={res_min_um:.3f}um  "
      f"taper_max={taper_max*1e6:.0f}um", flush=True)

if AMR:
    cap.enable_mesh_refinement(
        num_iterations=AMR_ITERS, max_DoFs=AMR_MAX_DOFS, tolerance=AMR_TOL,
        Dorfler_marking_fraction=AMR_DORFLER,
        save_iterations_data=True, save_iterations_mesh=False, nonconformal=False)
    print(f"AMR enabled: iters={AMR_ITERS} tol={AMR_TOL} max_DoFs={AMR_MAX_DOFS:,} "
          f"Dorfler={AMR_DORFLER}", flush=True)

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

import glob
import pandas as pd
import numpy as np


def _read_C_tg(csv_path):
    """trace<->ground off-diagonal capacitance (F) from a terminal-C.csv."""
    C = pd.read_csv(csv_path).to_numpy()[:, 1:].astype(float)
    gnd = int(np.argmax(np.diag(C)))
    trc = 1 - gnd
    return abs(C[trc, gnd]), C

# Per-AMR-iteration convergence of C itself (not just Palace's internal FEM
# error estimator) -- same pattern as squid_array_coupling/run_amr.py's
# per-iteration resonator-frequency tracking.
convergence = []
if AMR:
    for it_dir in sorted(glob.glob(os.path.join(cap._output_data_dir, "iteration*"))):
        csv_it = os.path.join(it_dir, "terminal-C.csv")
        if os.path.exists(csv_it):
            try:
                c_it, _ = _read_C_tg(csv_it)
                convergence.append({"iteration": os.path.basename(it_dir), "C_fF": c_it * 1e15})
            except Exception as e:
                convergence.append({"iteration": os.path.basename(it_dir), "error": str(e)})
    for i in range(1, len(convergence)):
        if "C_fF" in convergence[i] and "C_fF" in convergence[i - 1]:
            convergence[i]["rel_change_vs_prev"] = (
                abs(convergence[i]["C_fF"] - convergence[i - 1]["C_fF"]) / convergence[i]["C_fF"])
    print("AMR convergence (C_fF per iteration): "
          + json.dumps(convergence, indent=2), flush=True)

csv_path = cap._output_data_dir + "/terminal-C.csv"
C_tg, C = _read_C_tg(csv_path)
Cprime_F_per_m = C_tg / (LENGTH_UM * 1e-6)          # F/m
Cprime_fF_per_um = Cprime_F_per_m * 1e9              # 1 F/m = 1e9 fF/um

try:
    sizes = cap.retrieve_simulation_sizes()
except Exception:
    sizes = {}

row = {
    "tag": TAG, "dielectric_material": DIELECTRIC_MATERIAL,
    "cpw_width_um": CPW_WIDTH_UM, "cpw_gap_um": CPW_GAP_UM,
    "length_um": LENGTH_UM, "full3d": FULL3D,
    "metal_thickness_nm": METAL_THICKNESS_NM if FULL3D else None,
    "C_total_F": C_tg, "C_total_fF": C_tg * 1e15,
    "Cprime_F_per_m": Cprime_F_per_m, "Cprime_fF_per_um": Cprime_fF_per_um,
    "maxwell_matrix_F": C.tolist(),
    "dof": sizes.get("DoF"), "mesh_elements": sizes.get("MeshElements"),
    "amr": {"enabled": AMR, "iterations_requested": AMR_ITERS, "tolerance": AMR_TOL,
           "max_dofs": AMR_MAX_DOFS, "convergence": convergence} if AMR else {"enabled": False},
    "solve_s": solve_s, "wall_s": time.time() - t0,
}
with open(f"{OUTDIR}/results/{TAG}.json", "w") as f:
    json.dump(row, f, indent=2)

print(json.dumps({k: v for k, v in row.items() if k != "maxwell_matrix_F"}, indent=2),
      flush=True)

for pat in (f"{cap._output_data_dir}/../*.msh", f"{cap._output_data_dir}/paraview"):
    for p in glob.glob(pat):
        os.system(f"rm -rf '{p}'")
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
