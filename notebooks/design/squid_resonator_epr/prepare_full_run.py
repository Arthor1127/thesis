"""
Build + mesh + write the Palace config for the "full calculation" (finer
mesh, matching what was set in run_coarse_test.py), WITHOUT calling
eigen_sim.run(). Lets us inspect mesh size before committing to Palace's
(potentially slow) eigensolve. Run with `python -u` so prints aren't lost
if the process crashes hard (GMSH/Qt can SIGABRT without a catchable
Python exception).
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import (build_design, L_TRANSMON_NH, L_PROBE_NH, RES_LENGTH_MM,
                           CJ_TRANSMON_FF, CJ_PROBE_FF, QUBIT_COMPONENT_NAME,
                           JJ_INDEX_SINGLE, JJ_INDEX_SQUID_A)

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "runs/full_run1_prep"
PALACE_BIN = os.environ["PALACE_BIN"]
NUM_CPUS = int(os.environ.get("NSLOTS", 16))
os.makedirs(OUTDIR, exist_ok=True)

t0 = time.time()
print("[1/5] Building design...", flush=True)
design = build_design()
print(f"      OK ({time.time()-t0:.1f}s). Components: {list(design.components.keys())}", flush=True)

print("[2/5] design.rebuild()...", flush=True)
t1 = time.time()
design.rebuild()
print(f"      OK ({time.time()-t1:.1f}s)", flush=True)

print("[3/5] Importing SQDMetal PALACE eigenmode module...", flush=True)
from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.Utilities.Materials import MaterialInterface
import matplotlib
matplotlib.use("Agg", force=True)
print("      OK", flush=True)

print("[4/5] Setting up eigenmode simulation object + mesh (production settings, "
      "finer meshing from run_coarse_test.py)...", flush=True)
t2 = time.time()
user_defined_options = {
    "mesh_refinement": 0,
    "dielectric_material": "silicon",
    "starting_freq": 7.5e9,
    "number_of_freqs": 6,
    "solns_to_save": 6,
    "solver_order": 2,
    "solver_tol": 1.0e-8,
    "solver_maxits": 200,
    "fillet_resolution": 12,
    "palace_dir": PALACE_BIN,
    "num_cpus": NUM_CPUS,
}
eigen_sim = PALACE_Eigenmode_Simulation(
    name="full_run1_prep",
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

eigen_sim.fine_mesh_components(
    ["readout_1"],
    min_size=4e-6, max_size=150e-6, taper_dist_min=10e-6, metals_only=False,
)

# Q1_SQUID split in two, instead of one fine_mesh_components() pass over
# the whole component: the SQUID branches (l1+l2=600um long feed+tail
# leads either side of the junctions) make the component's total
# footprint ~611x37.5um, so a single min_size=2.5um pass was forcing that
# resolution across the whole strip AND the charge-island pad, when the
# pad's only real feature is pad_gap=10um (needs ~5um, not 2.5um) - see
# 2026-08-24 chat: this was the dominant cost in the mesh-generation OOM
# on the local workstation.
#
# Charge-island pad: circle of pad_radius=90um centered at the qubit's
# local origin (pos_x=pos_y=0, orientation=0, so local == global here).
# min_size=5um = ~half of pad_gap=10um, matching this repo's existing
# "half the smallest gap resolves it" convention (see
# --qubit-bg-min-size-um in build_and_eigenmode.py).
PAD_MARGIN_UM = 20e-6
PAD_RADIUS_UM = 90e-6
eigen_sim.fine_mesh_in_rectangle(
    -(PAD_RADIUS_UM + PAD_MARGIN_UM), -(PAD_RADIUS_UM + PAD_MARGIN_UM),
    (PAD_RADIUS_UM + PAD_MARGIN_UM), (PAD_RADIUS_UM + PAD_MARGIN_UM),
    min_size=5e-6, max_size=150e-6, taper_dist_min=5e-6, taper_dist_max=20e-6,
)

# SQUID branches: chord endpoints at pad_radius=90um, theta_2=180deg +-
# delta_angle=7.5deg, branches run in -x out to l1+l2+2*g1=608um past
# the chord, +-(w1/2+g1)=7um either side of each branch centerline.
# Computed bbox (see chat): x=[-697.23,-86.23]um, y=[-18.75,18.75]um;
# padded a little below. min_size kept at 2.5um - the user confirmed
# this resolution is adequate for the SQUID arms themselves.
eigen_sim.fine_mesh_in_rectangle(
    -710e-6, -25e-6, -80e-6, 25e-6,
    min_size=2.5e-6, max_size=50e-6, taper_dist_min=5e-6, taper_dist_max=20e-6,
)

for port in eigen_sim._ports[-2:]:
    coords = port['portCoords']
    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    eigen_sim.fine_mesh_in_rectangle(
        min(xs), min(ys), max(xs), max(ys),
        min_size=1.5e-6, max_size=5.0e-6,
        taper_dist_min=5e-6, taper_dist_max=20e-6,
    )
eigen_sim.setup_EPR_interfaces(
    metal_air=MaterialInterface("Aluminium-Vacuum"),
    substrate_air=MaterialInterface("Silicon-Vacuum"),
    substrate_metal=MaterialInterface("Silicon-Aluminium"),
)
print(f"      OK ({time.time()-t2:.1f}s)", flush=True)

print("[5/5] prepare_simulation() - writing mesh + Palace config (no solve)...", flush=True)
t3 = time.time()
eigen_sim.prepare_simulation()
print(f"      OK ({time.time()-t3:.1f}s)", flush=True)
print(f"\n=== PREP COMPLETE in {time.time()-t0:.1f}s. Files written to {OUTDIR} ===", flush=True)
