"""
Electrostatic capacitance of a SINGLE circular metal island (radius a) sitting
in a ground plane with an annular vacuum gap of width g (outer radius b = a + g).
No junction, no coupler -- just island <-> ground.

Full-3D: the metal is given a real vertical thickness t so the island-edge /
ground-edge "coaxial cylinder" capacitance is physically present. This is the
sweep needed to test

    C ?= 2*pi * eps_eff * t / ln(b/a),   b = a + g,   eps_eff = eps0*(1+eps_r)/2

Env: RADIUS_UM, GAP_UM, THICK_NM, [EPS_R=11.9], PALACE_BIN, [NSLOTS=4],
     [OUTDIR=./island_cap], [TAG]
Writes: {OUTDIR}/results/{TAG}.json   (one row of the sweep)
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

RADIUS_UM = float(os.environ["RADIUS_UM"])
GAP_UM    = float(os.environ["GAP_UM"])
THICK_NM  = float(os.environ["THICK_NM"])
EPS_R     = float(os.environ.get("EPS_R", "11.9"))
OUTDIR    = os.environ.get("OUTDIR", "./island_cap")
PALACE_BIN = os.environ["PALACE_BIN"]
NCPU      = int(os.environ.get("NSLOTS", "4"))
TAG       = os.environ.get("TAG", f"r{RADIUS_UM:g}_g{GAP_UM:g}_t{THICK_NM:g}")
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

print(f"[{TAG}] radius={RADIUS_UM}um gap={GAP_UM}um thickness={THICK_NM}nm "
      f"eps_r={EPS_R} ncpu={NCPU}", flush=True)

t0 = time.time()
from qiskit_metal import designs, draw, Dict
from qiskit_metal.qlibrary.core import QComponent
# qiskit_metal flips the mpl backend to QtAgg on import; force it back to Agg
# AFTER that import and BEFORE the SQDMetal PALACE import, or SQDMetal's headless
# plot calls (display_conductor_indices) crash with a QtAgg ImportError.
import matplotlib
matplotlib.use("Agg", force=True)


class CircularIsland(QComponent):
    """A bare circular metal disc (radius `radius`) with a concentric circular
    ground-plane cutout (radius `radius` + `gap`). Emits exactly two polys:
    the disc (metal) and the pocket (subtract) -- nothing else."""
    default_options = Dict(pos_x="0um", pos_y="0um", chip="main",
                           radius="90um", gap="10um", resolution="256")

    def make(self):
        p = self.p
        res = int(p.resolution)
        disc = draw.Point(0, 0).buffer(p.radius, resolution=res)
        pocket = draw.Point(0, 0).buffer(p.radius + p.gap, resolution=res)
        disc, pocket = draw.translate([disc, pocket], p.pos_x, p.pos_y)
        self.add_qgeometry("poly", {"island": disc}, chip=p.chip)
        self.add_qgeometry("poly", {"pocket": pocket}, subtract=True, chip=p.chip)


design = designs.DesignPlanar({}, overwrite_enabled=True)
# Chip large enough that the outer PEC wall does not perturb the island->ground
# capacitance: half-width = b + 4*a  (>= 4 island radii past the gap -- the
# field past the gap decays fast, boundary shift is <0.5%). A 10*a margin (used
# earlier) bloated the domain volume ~6x for no accuracy gain.
half_um = (RADIUS_UM + GAP_UM) + 4.0 * RADIUS_UM
design.chips.main.size.size_x = f"{2*half_um/1000:.5f}mm"
design.chips.main.size.size_y = f"{2*half_um/1000:.5f}mm"
design.chips.main.size.size_z = "500um"
design.chips.main.size.center_x = "0mm"
design.chips.main.size.center_y = "0mm"

CircularIsland(design, "isl", options=dict(
    radius=f"{RADIUS_UM}um", gap=f"{GAP_UM}um", resolution=200))
design.rebuild()
print(f"geometry built  ({time.time()-t0:.1f}s)  chip={design.chips.main.size.size_x}", flush=True)

from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation

BG_MAX_UM = float(os.environ.get("BG_MAX_UM", "45"))
cap = PALACE_Capacitance_Simulation(
    name=f"ci_{TAG}", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "dielectric_material": "silicon",
        # Order 1 is well-converged for a capacitance extraction; order 2 gave
        # 5-8x the DoF (and 25-68 GB / task, some blowing to 27M unknowns) for
        # a sub-percent change in C. Do NOT raise this without a memory guard.
        "solver_order": int(os.environ.get("SOLVER_ORDER", "1")),
        "solver_tol": 1e-8,
        "solver_maxits": 200,
        "solns_to_save": 0,
        "fillet_resolution": 10,
        "mesh_max": BG_MAX_UM * 1e-6,     # cap the bulk mesh
        "palace_dir": PALACE_BIN,
        "num_cpus": NCPU,
    })
cap.add_metallic(1)
cap.add_ground_plane()

# FULL3D=1 -> metal has a real vertical thickness t (needed for any genuine
# t-dependence, but GMSH must put >=1 element across t, so min_size ~ t; at
# t<~0.2um on ~100um features this explodes / hangs). Default: quasi-3D (metal
# as a zero-thickness sheet) -- fast and exact for the in-plane capacitance;
# the coax-edge term (the only t-dependent piece) is then simply absent, which
# is itself the result: it is a <1% correction here.
FULL3D = os.environ.get("FULL3D", "0") == "1"
if FULL3D:
    cap.enforce_full_3D_simulation(metal_thickness=THICK_NM * 1e-9, substrate_trenching=0.0)
    res_min = max(THICK_NM * 1e-9, 0.15e-6)          # >= metal thickness
elif os.environ.get("RES_MIN_UM"):
    res_min = float(os.environ["RES_MIN_UM"]) * 1e-6          # explicit override
else:
    # ~4 elements across the gap; floor 1.0 um. Going below ~1 um on the thin
    # curved ring makes GMSH thrash on sliver elements (g=3 um at 0.6 um hung
    # in meshing for >10 min). ~4 elements across the gap should give C to a
    # few % at order 1; a res_min convergence check on g=3 is worthwhile.
    res_min = max(GAP_UM / 4.0, 1.0) * 1e-6

# fine mesh a thin ring around the rim only.
taper_max = max(3.0 * GAP_UM, 14.0) * 1e-6
cap.fine_mesh_components(["isl"], min_size=res_min, max_size=BG_MAX_UM * 1e-6,
                         taper_dist_min=1.5 * GAP_UM * 1e-6, taper_dist_max=taper_max,
                         metals_only=False)
print(f"mesh: FULL3D={FULL3D}  res_min={res_min*1e6:.3f}um  taper_max={taper_max*1e6:.0f}um", flush=True)

tp = time.time()
cap.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

ts = time.time()
try:
    cap.run()                          # also runs retrieve_data() -> may hit a
except Exception as e:                  # headless-mpl plot error AFTER the solve
    print(f"cap.run() post-processing raised (solve output still written): {e}", flush=True)
solve_s = time.time() - ts
print(f"SOLVE: {solve_s:.1f} s", flush=True)

# read the Maxwell matrix straight from Palace's CSV -- independent of the
# plotting step that may have crashed.
import pandas as pd
csv_path = cap._output_data_dir + "/terminal-C.csv"
C = pd.read_csv(csv_path).to_numpy()[:, 1:].astype(float)   # strip index column
# 2 conductors: island + ground plane. island<->ground cap = |off-diagonal|.
# (identify ground = larger self-capacitance, just for the record)
gnd = int(np.argmax(np.diag(C)))
isl = 1 - gnd
C_ig = abs(C[isl, gnd])               # F
C_isl_self = abs(C[isl, isl])

# analytic "coaxial edge" prediction
eps0 = 8.8541878128e-12
eps_eff = eps0 * (1.0 + EPS_R) / 2.0
a = RADIUS_UM * 1e-6
b = (RADIUS_UM + GAP_UM) * 1e-6
t = THICK_NM * 1e-9
C_formula = 2.0 * np.pi * eps_eff * t / np.log(b / a)

try:
    sizes = cap.retrieve_simulation_sizes()
except Exception:
    sizes = {}

row = {
    "tag": TAG,
    "radius_um": RADIUS_UM, "gap_um": GAP_UM, "thickness_nm": THICK_NM, "eps_r": EPS_R,
    "b_over_a": b / a,
    "C_island_ground_F": C_ig,
    "C_island_ground_fF": C_ig * 1e15,
    "C_island_self_fF": C_isl_self * 1e15,
    "C_formula_coax_fF": C_formula * 1e15,
    "ratio_sim_over_formula": C_ig / C_formula if C_formula > 0 else None,
    "maxwell_matrix_F": C.tolist(),
    "dof": sizes.get("DoF"), "mesh_elements": sizes.get("MeshElements"),
    "solve_s": solve_s, "wall_s": time.time() - t0,
}
with open(f"{OUTDIR}/results/{TAG}.json", "w") as f:
    json.dump(row, f, indent=2)

print(json.dumps({k: v for k, v in row.items() if k != "maxwell_matrix_F"}, indent=2), flush=True)
print(f"\nC(island->ground) = {C_ig*1e15:.4f} fF   |   coax formula = "
      f"{C_formula*1e15:.4f} fF   |   sim/formula = {C_ig/C_formula:.2f}", flush=True)

# drop the heavy per-run artifacts (mesh + paraview dir); keep only the CSVs + json
import glob
for pat in (f"{cap._output_data_dir}/../*.msh", f"{cap._output_data_dir}/paraview"):
    for p in glob.glob(pat):
        os.system(f"rm -rf '{p}'")
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
