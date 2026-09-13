"""
Capacitance of a 2-CircTransmon + 1 CapNInterdigitalRound pair, for the
coupler-geometry sweep. Extracts the direct island-island mutual
|C_Maxwell[Q1][Q2]| = the coupling capacitance C_c, plus each island's
capacitance to the ground plane.

WITH_COUPLER=0 draws only the two pads at the same separation -> the bare
pad-to-pad reference to subtract.

Palace AMR is used with a low max_DoFs cap (single coupler stays well under the
~9 M-DOF size where Palace's flux-recovery error estimator stalls).

Env: FINGER_COUNT, FINGER_LENGTH_UM, FINGER_GAP_UM, FINGER_WIDTH_UM,
     FINGER_TIP_GAP_UM, SPINE_WIDTH_UM, NECK_UM, PAD_RADIUS_UM, PAD_GAP_UM,
     WITH_COUPLER(1), AMR(1), AMR_ITERS(4), AMR_MAX_DOFS(4000000), AMR_TOL(5e-3),
     COUPLER_RES_UM(2.0), TRANSMON_RES_UM(14), SOLVER_ORDER(2), BG_MAX_UM(80),
     PALACE_BIN, NSLOTS, OUTDIR, TAG
No ParaView output.
"""
import os
import sys
import json
import time
import glob

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

FC   = int(os.environ.get("FINGER_COUNT", "8"))
FL   = float(os.environ.get("FINGER_LENGTH_UM", "80"))
FG   = float(os.environ.get("FINGER_GAP_UM", "4"))
FW   = float(os.environ.get("FINGER_WIDTH_UM", "8"))
FTG  = float(os.environ.get("FINGER_TIP_GAP_UM", "6"))
SW   = float(os.environ.get("SPINE_WIDTH_UM", "12"))
NECK = float(os.environ.get("NECK_UM", "60"))
PADR = float(os.environ.get("PAD_RADIUS_UM", "136"))
PADG = float(os.environ.get("PAD_GAP_UM", "20"))
WITH_COUPLER = os.environ.get("WITH_COUPLER", "1") == "1"

USE_AMR    = os.environ.get("AMR", "1") == "1"
AMR_ITERS  = int(os.environ.get("AMR_ITERS", "6"))
AMR_MAXD   = int(os.environ.get("AMR_MAX_DOFS", "1500000"))
AMR_TOL    = float(os.environ.get("AMR_TOL", "1e-3"))
AMR_DORFLER = float(os.environ.get("AMR_DORFLER", "0.35"))
CPL_RES   = float(os.environ.get("COUPLER_RES_UM", "2.0"))
TRM_RES   = float(os.environ.get("TRANSMON_RES_UM", "14"))
ORDER     = int(os.environ.get("SOLVER_ORDER", "2"))
BG_MAX    = float(os.environ.get("BG_MAX_UM", "80"))
OUTDIR    = os.environ.get("OUTDIR", "./coupler_sweep")
PALACE_BIN = os.environ["PALACE_BIN"]
NCPU      = int(os.environ.get("NSLOTS", "8"))
TAG       = os.environ.get("TAG", f"N{FC}_L{FL:g}_g{FG:g}_w{FW:g}_t{FTG:g}"
                           + ("" if WITH_COUPLER else "_nocpl"))
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

print(f"[{TAG}] FC={FC} FL={FL} FG={FG} FW={FW} FTG={FTG} neck={NECK} "
      f"with_coupler={WITH_COUPLER} AMR={USE_AMR} ncpu={NCPU}", flush=True)

t0 = time.time()
from qiskit_metal import designs, draw, Dict           # noqa: F401
from qiskit_metal.qlibrary.core import QComponent       # noqa: F401
import matplotlib
matplotlib.use("Agg", force=True)

import build_coupler_pair as bcp

design = bcp.build_coupler_pair(
    finger_count=FC, finger_length_um=FL, finger_gap_um=FG, finger_width_um=FW,
    finger_tip_gap_um=FTG, spine_width_um=SW, neck_um=NECK,
    pad_radius_um=PADR, pad_gap_um=PADG, with_coupler=WITH_COUPLER)
design.rebuild()
comps = list(design.components.keys())
transmons = [c for c in comps if c.startswith("Q")]
couplers = [c for c in comps if c.startswith("coupler")]
print(f"geometry built ({time.time()-t0:.1f}s)  components={comps}  "
      f"chip={design.chips.main.size.size_x} x {design.chips.main.size.size_y}", flush=True)

from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation

cap = PALACE_Capacitance_Simulation(
    name=f"cp_{TAG}", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "mesh_refinement": 0,
        "dielectric_material": "silicon",
        "solver_order": ORDER,
        "solver_tol": 1e-8,
        "solver_maxits": 300,
        "solns_to_save": 0,
        "fillet_resolution": 12,
        "mesh_max": BG_MAX * 1e-6,
        "palace_dir": PALACE_BIN,
        "num_cpus": NCPU,
    })
cap.add_metallic(1)
cap.add_ground_plane()
cap.fine_mesh_components(transmons, min_size=TRM_RES * 1e-6, max_size=BG_MAX * 1e-6,
                         taper_dist_min=12e-6, taper_dist_max=60e-6, metals_only=False)
if couplers:
    cap.fine_mesh_components(couplers, min_size=CPL_RES * 1e-6, max_size=BG_MAX * 1e-6,
                             taper_dist_min=5e-6, taper_dist_max=30e-6, metals_only=False)
if USE_AMR:
    cap.enable_mesh_refinement(num_iterations=AMR_ITERS, max_DoFs=AMR_MAXD, tolerance=AMR_TOL,
                               Dorfler_marking_fraction=AMR_DORFLER,
                               save_iterations_data=True, save_iterations_mesh=False,
                               nonconformal=False)

tp = time.time()
cap.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

ts = time.time()
try:
    cap.run()
except Exception as e:
    print(f"cap.run() post-processing raised: {e}", flush=True)
solve_s = time.time() - ts
print(f"SOLVE: {solve_s:.1f} s", flush=True)

import pandas as pd
odir = cap._output_data_dir


def read_C(path):
    raw = pd.read_csv(path)
    return [c.strip() for c in raw.columns[1:]], raw.to_numpy()[:, 1:].astype(float)


def summarise(names, C):
    gnd = int(np.argmax(np.diag(C)))
    n = C.shape[0]
    isl = [i for i in range(n) if i != gnd]
    out = {"n_conductors": n, "terminal_names": names, "ground_index": gnd,
           "C_mutual_all_fF": {f"{names[a]}--{names[b]}": abs(C[a, b]) * 1e15
                               for a in range(n) for b in range(a + 1, n)},
           "C_to_ground_fF": {names[i]: abs(C[i, gnd]) * 1e15 for i in isl},
           "C_selfdiag_fF": {names[i]: C[i, i] * 1e15 for i in isl}}
    # the coupling capacitance = largest island-island mutual
    mm = [(abs(C[a, b]) * 1e15, a, b) for a in isl for b in isl if a < b]
    if mm:
        out["C_c_fF"] = max(mm)[0]
    return out


names, C = read_C(f"{odir}/terminal-C.csv")
final = summarise(names, C)

cycles = []
for cdir in sorted(glob.glob(f"{odir}/**/terminal-C.csv", recursive=True)):
    if os.path.abspath(cdir) == os.path.abspath(f"{odir}/terminal-C.csv"):
        continue
    try:
        cn, cc = read_C(cdir)
        s = summarise(cn, cc)
        cycles.append({"iter": os.path.basename(os.path.dirname(cdir)),
                       "n_conductors": s["n_conductors"], "C_c_fF": s.get("C_c_fF")})
    except Exception:
        pass

try:
    sizes = cap.retrieve_simulation_sizes()
except Exception:
    sizes = {}

row = {
    "tag": TAG, "with_coupler": WITH_COUPLER,
    "finger_count": FC, "finger_length_um": FL, "finger_gap_um": FG,
    "finger_width_um": FW, "finger_tip_gap_um": FTG, "spine_width_um": SW,
    "neck_um": NECK, "pad_radius_um": PADR, "pad_gap_um": PADG,
    "eps_r_sim": 11.45, "solver_order": ORDER, "amr": USE_AMR,
    **final,
    "maxwell_C_matrix_F": C.tolist(),
    "amr_cycles": cycles,
    "dof": sizes.get("DoF"), "mesh_elements": sizes.get("MeshElements"),
    "solve_s": solve_s, "wall_s": time.time() - t0,
}
with open(f"{OUTDIR}/results/{TAG}.json", "w") as f:
    json.dump(row, f, indent=2)
for src in glob.glob(f"{odir}/terminal-*.csv"):
    os.system(f"cp '{src}' '{OUTDIR}/results/{TAG}_{os.path.basename(src)}'")

print(json.dumps({k: v for k, v in row.items() if k != "maxwell_C_matrix_F"},
                 indent=2, default=str)[:4000], flush=True)
for pat in (f"{odir}/../*.msh", f"{odir}/paraview"):
    for p in glob.glob(pat):
        os.system(f"rm -rf '{p}'")
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
