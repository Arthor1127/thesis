"""
Adaptive mesh refinement (AMR) run: start from a modest, UNBIASED uniform
mesh (single fine_mesh_components pass over readout_1 + Q1_SQUID as a
whole, no manual pad/squid_arm/junction split - that was our own guess;
this time we let Palace's own error indicator decide where refinement
actually matters), then let Palace's AMR (SaveAdaptMesh=True,
SaveAdaptIterations=True) refine over a few iterations.

Unlike the earlier probe scripts, this calls eigen_sim.run() and lets it
actually solve (AMR needs a real solve + error estimate each iteration to
decide where to refine) - so this is NOT a kill-early probe. Wrapped in a
bash `timeout` at the SGE-job level since nobody will be watching this
session live (see 2026-08-25 chat).

After it finishes (or times out), inspect outputFiles/iteration1/,
iteration2/, ... (per-AMR-iteration ParaView output) for the per-element
`Indicator` field - that's Palace's own error indicator, i.e. exactly
where it decided to refine. See docs/src/guide/postprocessing.md in the
Palace source (~/repo/palace-src) for the field reference.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import (build_design, L_TRANSMON_NH, L_PROBE_NH,
                           CJ_TRANSMON_FF, CJ_PROBE_FF, QUBIT_COMPONENT_NAME,
                           JJ_INDEX_SINGLE, JJ_INDEX_SQUID_A)

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "runs/amr_run1"
PALACE_BIN = os.environ["PALACE_BIN"]
NUM_CPUS = int(os.environ.get("NSLOTS", 16))
AMR_ITERS = int(os.environ.get("AMR_ITERS", 2))
os.makedirs(OUTDIR, exist_ok=True)

t0 = time.time()
print(f"[1/6] Building design...", flush=True)
design = build_design()
design.rebuild()
print(f"      OK ({time.time()-t0:.1f}s). Components: {list(design.components.keys())}", flush=True)

print("[2/6] Importing SQDMetal PALACE eigenmode module...", flush=True)
from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.Utilities.Materials import MaterialInterface
import matplotlib
matplotlib.use("Agg", force=True)
print("      OK", flush=True)

print("[3/6] Setting up eigenmode simulation object...", flush=True)
t1 = time.time()
user_defined_options = {
    "mesh_refinement": 0,   # UniformLevels; AMR set separately below via enable_mesh_refinement()
    "dielectric_material": "silicon",
    "starting_freq": 7.0e9,
    "number_of_freqs": 2,
    "solns_to_save": 2,
    "solver_order": 1,
    "solver_tol": 1.0e-5,
    "solver_maxits": 300,
    "fillet_resolution": 10,
    "palace_dir": PALACE_BIN,
    "num_cpus": NUM_CPUS,
}
eigen_sim = PALACE_Eigenmode_Simulation(
    name="amr_run1",
    metal_design=design,
    sim_parent_directory=OUTDIR,
    mode="simPC",
    meshing="GMSH",
    user_options=user_defined_options,
    create_files=True,
)
eigen_sim.add_metallic(1)
eigen_sim.add_ground_plane()
eigen_sim.create_port_JosephsonJunction(
    QUBIT_COMPONENT_NAME, junction_index=JJ_INDEX_SINGLE,
    L_J=L_TRANSMON_NH * 1e-9, C_J=CJ_TRANSMON_FF * 1e-15)
eigen_sim.create_port_JosephsonJunction(
    QUBIT_COMPONENT_NAME, junction_index=JJ_INDEX_SQUID_A,
    L_J=L_PROBE_NH * 1e-9, C_J=CJ_PROBE_FF * 1e-15)

# Taper sizes updated per the 2026-08-25 AMR-1 result: the Indicator
# field showed refinement concentrated around the readout resonator
# (top-indicator points along its trace edges, over a wide span) and
# barely touching the qubit pad/SQUID arms (~1-4% of indicator mass).
# So the readout's fine-mesh zone is widened and the qubit's narrowed to
# match - same direction as the run_coarse_test.py update, applied here
# to the base mesh AMR then refines further from.
print("[4/6] Base fine mesh (readout widened, qubit narrowed per AMR-1 findings)...", flush=True)
t2 = time.time()
eigen_sim.fine_mesh_components(
    ["readout_1"],
    min_size=4e-6, max_size=150e-6, taper_dist_min=40e-6, metals_only=False,
)
eigen_sim.fine_mesh_components(
    [QUBIT_COMPONENT_NAME],
    min_size=2.5e-6, max_size=150e-6, taper_dist_min=8e-6, metals_only=False,
)
eigen_sim.setup_EPR_interfaces(
    metal_air=MaterialInterface("Aluminium-Vacuum"),
    substrate_air=MaterialInterface("Silicon-Vacuum"),
    substrate_metal=MaterialInterface("Silicon-Aluminium"),
)
print(f"      OK ({time.time()-t2:.1f}s)", flush=True)

print(f"[5/6] Enabling AMR: {AMR_ITERS} iterations, tol=1e-2, saving per-iteration "
      "mesh+data (Indicator field)...", flush=True)
eigen_sim.enable_mesh_refinement(
    num_iterations=AMR_ITERS, tolerance=1e-2, Dorfler_marking_fraction=0.7,
    save_iterations_data=True, save_iterations_mesh=True, nonconformal=False,
)

print("[6/6] prepare_simulation() + run() (full solve, AMR needs real error "
      "estimates each iteration - this is NOT a kill-early probe)...", flush=True)
t3 = time.time()
eigen_sim.prepare_simulation()
eigen_sim.run()
print(f"      OK ({time.time()-t3:.1f}s)", flush=True)

port_epr = eigen_sim.retrieve_mode_port_EPR()
mesh_sizes = eigen_sim.retrieve_simulation_sizes()
print(f"\n=== AMR RUN COMPLETE in {time.time()-t0:.1f}s ===", flush=True)
print(f"Final mesh: {mesh_sizes}", flush=True)
print(f"Eigenfrequencies (GHz): {[f/1e9 for f in port_epr['eigenfrequencies']]}", flush=True)
print(f"Loaded Q: {list(port_epr['loaded_Q'])}", flush=True)
print(f"mat_mode_port:\n{port_epr['mat_mode_port']}", flush=True)
