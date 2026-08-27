"""Mesh-only check matching build_and_eigenmode.py's current exact
default settings (readout taper widened, junction min raised to 6um,
squid-loop blanket rectangle removed) - to see element count before
committing to a full Palace solve."""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import (build_design, L_TRANSMON_NH, L_PROBE_NH,
                           CJ_TRANSMON_FF, CJ_PROBE_FF, QUBIT_COMPONENT_NAME,
                           JJ_INDEX_SINGLE, JJ_INDEX_SQUID_A)

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "runs/check_prod_mesh"
PALACE_BIN = os.environ.get("PALACE_BIN", "/nonexistent")
os.makedirs(OUTDIR, exist_ok=True)

t0 = time.time()
design = build_design()
design.rebuild()
print(f"[1/3] design OK ({time.time()-t0:.1f}s)", flush=True)

from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.Utilities.Materials import MaterialInterface
import matplotlib
matplotlib.use("Agg", force=True)

user_defined_options = {
    "mesh_refinement": 0, "dielectric_material": "silicon",
    "starting_freq": 7.0e9, "number_of_freqs": 3, "solns_to_save": 3,
    "solver_order": 1, "solver_tol": 1.0e-5, "solver_maxits": 300,
    "fillet_resolution": 12, "palace_dir": PALACE_BIN, "num_cpus": 1,
}
eigen_sim = PALACE_Eigenmode_Simulation(
    name="checkprod", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", user_options=user_defined_options, create_files=True,
)
eigen_sim.add_metallic(1)
eigen_sim.add_ground_plane()
eigen_sim.create_port_JosephsonJunction(
    QUBIT_COMPONENT_NAME, junction_index=JJ_INDEX_SINGLE,
    L_J=L_TRANSMON_NH * 1e-9, C_J=CJ_TRANSMON_FF * 1e-15)
eigen_sim.create_port_JosephsonJunction(
    QUBIT_COMPONENT_NAME, junction_index=JJ_INDEX_SQUID_A,
    L_J=L_PROBE_NH * 1e-9, C_J=CJ_PROBE_FF * 1e-15)

# Matches build_and_eigenmode.py defaults exactly (bg_min=6, bg_max=150,
# qubit_min=2.0, qubit_max=20, junction_min=6, junction_max=20)
eigen_sim.fine_mesh_components(
    ["readout_1"], min_size=6e-6, max_size=150e-6, taper_dist_min=40e-6, metals_only=False,
)
eigen_sim.fine_mesh_components(
    [QUBIT_COMPONENT_NAME], min_size=2.0e-6, max_size=150e-6, taper_dist_min=5e-6, metals_only=False,
)
# junction-local fine_mesh_in_rectangle removed - isolated as the cause
# of the mesh blowup (2026-08-25)
eigen_sim.setup_EPR_interfaces(
    metal_air=MaterialInterface("Aluminium-Vacuum"),
    substrate_air=MaterialInterface("Silicon-Vacuum"),
    substrate_metal=MaterialInterface("Silicon-Aluminium"),
)
print("[2/3] mesh setup OK, calling prepare_simulation()...", flush=True)
t1 = time.time()
eigen_sim.prepare_simulation()
print(f"[3/3] prepare_simulation() OK ({time.time()-t1:.1f}s)", flush=True)

msh_path = os.path.join(OUTDIR, "checkprod", "checkprod.msh")
if os.path.exists(msh_path):
    with open(msh_path) as f:
        for line in f:
            if line.strip() == "$Elements":
                print(f"MESH_ELEMENTS={int(next(f).strip())}", flush=True)
                break
print(f"\n=== DONE in {time.time()-t0:.1f}s ===", flush=True)
