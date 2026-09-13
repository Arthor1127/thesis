"""
Capacitance matrix of the 6-transmon interdigital-coupler chain with
Palace ADAPTIVE MESH REFINEMENT (AMR).

The manual coupler-mesh sweep (3.0 / 2.0 / 1.5 um) does not converge the
interdigital-coupler capacitances -- they fall ~13 % per refinement because the
metal-edge field singularity in the 4 um finger gaps is under-resolved, and
GMSH chokes below ~1 um on the rounded finger tips. AMR lets Palace add
elements exactly where its a-posteriori error estimate is largest (the finger
edges) and reports the error indicator per cycle, giving a real convergence
criterion.

Starts from a modest uniform coupler mesh (COUPLER_RES_UM, default 2.5) and
refines for AMR_ITERS cycles to tolerance AMR_TOL, capped at AMR_MAX_DOFS
(RAM guard). Per-cycle terminal-C.csv files are collected so the matrix can be
watched converging.

Env: N_TRANSMONS(6), COUPLER_RES_UM(2.5), TRANSMON_RES_UM(12), SOLVER_ORDER(2),
     BG_MAX_UM(80), MARGIN_X_MM(1.2), MARGIN_Y_MM(1.2),
     AMR_ITERS(6), AMR_TOL(5e-3), AMR_DORFLER(0.7), AMR_MAX_DOFS(12000000),
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

N          = int(os.environ.get("N_TRANSMONS", "6"))
CPL_RES    = float(os.environ.get("COUPLER_RES_UM", "2.5"))
TRM_RES    = float(os.environ.get("TRANSMON_RES_UM", "12.0"))
ORDER      = int(os.environ.get("SOLVER_ORDER", "2"))
BG_MAX_UM  = float(os.environ.get("BG_MAX_UM", "80"))
MARGIN_X   = float(os.environ.get("MARGIN_X_MM", "1.2"))
MARGIN_Y   = float(os.environ.get("MARGIN_Y_MM", "1.2"))
AMR_ITERS  = int(os.environ.get("AMR_ITERS", "6"))
AMR_TOL    = float(os.environ.get("AMR_TOL", "5e-3"))
AMR_DORFLER = float(os.environ.get("AMR_DORFLER", "0.7"))
AMR_MAX_DOFS = int(os.environ.get("AMR_MAX_DOFS", "12000000"))
OUTDIR     = os.environ.get("OUTDIR", "./chain_cap")
PALACE_BIN = os.environ["PALACE_BIN"]
NCPU       = int(os.environ.get("NSLOTS", "8"))
TAG        = os.environ.get("TAG", f"chain{N}_amr_o{ORDER}")
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

print(f"[{TAG}] N={N} base coupler_res={CPL_RES}um transmon_res={TRM_RES}um "
      f"order={ORDER} | AMR iters={AMR_ITERS} tol={AMR_TOL} max_DoFs={AMR_MAX_DOFS:,} "
      f"ncpu={NCPU}", flush=True)

t0 = time.time()
from qiskit_metal import designs, draw, Dict            # noqa: F401
from qiskit_metal.qlibrary.core import QComponent        # noqa: F401
import matplotlib
matplotlib.use("Agg", force=True)

import build_chain
from build_design import _design_bounds_mm

design = build_chain.build_chain(n_transmons=N)
minx, miny, maxx, maxy = _design_bounds_mm(design)
sx, sy = (maxx - minx) + MARGIN_X, (maxy - miny) + MARGIN_Y
design.chips.main.size.size_x = f"{sx}mm"
design.chips.main.size.size_y = f"{sy}mm"
design.chips.main.size.center_x = f"{0.5*(minx+maxx)}mm"
design.chips.main.size.center_y = f"{0.5*(miny+maxy)}mm"
design.rebuild()

comps = list(design.components.keys())
transmons = [c for c in comps if c.startswith("Q")]
couplers  = [c for c in comps if c.startswith("coupler")]
print(f"geometry built ({time.time()-t0:.1f}s)  chip {sx:.3f} x {sy:.3f} mm", flush=True)

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
        "solns_to_save": 0,
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
cap.enable_mesh_refinement(
    num_iterations=AMR_ITERS, max_DoFs=AMR_MAX_DOFS, tolerance=AMR_TOL,
    Dorfler_marking_fraction=AMR_DORFLER,
    save_iterations_data=True, save_iterations_mesh=False, nonconformal=False)

tp = time.time()
cap.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

ts = time.time()
try:
    cap.run()
except Exception as e:
    print(f"cap.run() post-processing raised (solve output still written): {e}", flush=True)
solve_s = time.time() - ts
print(f"SOLVE (all AMR cycles): {solve_s:.1f} s", flush=True)

import pandas as pd
odir = cap._output_data_dir


def read_C(path):
    raw = pd.read_csv(path)
    names = [c.strip() for c in raw.columns[1:]]
    return names, raw.to_numpy()[:, 1:].astype(float)


def matrix_summary(names, C):
    gnd = int(np.argmax(np.diag(C)))
    n = C.shape[0]
    isl = [i for i in range(n) if i != gnd]
    return {
        "ground_index": gnd,
        "C_island_to_groundplane_fF": {names[i]: abs(C[i, gnd]) * 1e15 for i in isl},
        "C_island_selfdiag_fF": {names[i]: C[i, i] * 1e15 for i in isl},
        "C_mutual_fF": {f"{names[a]}--{names[b]}": abs(C[a, b]) * 1e15
                        for a in range(n) for b in range(a + 1, n)},
    }


# final matrix
names, C = read_C(f"{odir}/terminal-C.csv")
final = matrix_summary(names, C)

# per-AMR-cycle matrices (Palace writes them under the iterations dir)
cycles = []
for cdir in sorted(glob.glob(f"{odir}/**/terminal-C.csv", recursive=True)):
    if os.path.abspath(cdir) == os.path.abspath(f"{odir}/terminal-C.csv"):
        continue
    try:
        cn, cc = read_C(cdir)
        cycles.append({"path": os.path.relpath(cdir, odir),
                       "NN_mutual_fF": [abs(cc[i, i + 1]) * 1e15 for i in range(1, cc.shape[0] - 1)],
                       "summary": matrix_summary(cn, cc)})
    except Exception:
        pass

# AMR error-indicator history
err_hist = None
for ename in ("error-indicators.csv", "domain-error-estimates.csv"):
    p = f"{odir}/{ename}"
    if os.path.exists(p):
        try:
            err_hist = pd.read_csv(p).to_dict(orient="list")
        except Exception:
            pass

try:
    sizes = cap.retrieve_simulation_sizes()
except Exception:
    sizes = {}

row = {
    "tag": TAG, "n_transmons": N, "method": "palace_AMR",
    "base_coupler_res_um": CPL_RES, "transmon_res_um": TRM_RES,
    "solver_order": ORDER, "amr": {"iters": AMR_ITERS, "tol": AMR_TOL,
                                   "dorfler": AMR_DORFLER, "max_dofs": AMR_MAX_DOFS},
    "chip_mm": [sx, sy], "eps_r_sim": 11.45,
    "terminal_names": names,
    "final": final,
    "maxwell_C_matrix_F": C.tolist(),
    "amr_cycles": cycles,
    "amr_error_history": err_hist,
    "dof": sizes.get("DoF"), "mesh_elements": sizes.get("MeshElements"),
    "solve_s": solve_s, "wall_s": time.time() - t0,
}
with open(f"{OUTDIR}/results/{TAG}.json", "w") as f:
    json.dump(row, f, indent=2)
for src in glob.glob(f"{odir}/terminal-*.csv") + glob.glob(f"{odir}/error-indicators.csv"):
    os.system(f"cp '{src}' '{OUTDIR}/results/{TAG}_{os.path.basename(src)}'")

print(json.dumps({k: v for k, v in row.items()
                  if k not in ("maxwell_C_matrix_F",)}, indent=2, default=str)[:6000], flush=True)

for pat in (f"{odir}/../*.msh", f"{odir}/paraview"):
    for p in glob.glob(pat):
        os.system(f"rm -rf '{p}'")
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
