"""
Very coarse smoke test of the full build_design -> Palace eigenmode -> EPR
retrieval pipeline. Not physically meaningful (mesh is deliberately coarse,
few modes, no AMR) - the goal is only to confirm the pipeline runs end to
end and produces a well-formed epr_result.json.
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import build_design, _design_bounds_mm, QUBIT_COMPONENT_NAME, \
    JJ_INDEX_SINGLE, JJ_INDEX_SQUID_A, L_TRANSMON_NH, L_PROBE_NH, CJ_TRANSMON_FF, CJ_PROBE_FF


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--palace-bin", required=True)
    ap.add_argument("--num-cpus", type=int, default=16)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    t0 = time.time()
    print("[1/3] Building design (auto-sized chip from geometry bounds)...")
    design = build_design()
    minx, miny, maxx, maxy = _design_bounds_mm(design)
    print(f"      Design bounds (mm): X=[{minx:.4f},{maxx:.4f}]  Y=[{miny:.4f},{maxy:.4f}]")
    print(f"      Chip size (mm): {design.chips.main.size.size_x} x {design.chips.main.size.size_y}, "
          f"centered at ({design.chips.main.size.center_x}, {design.chips.main.size.center_y})")
    print(f"      Components: {list(design.components.keys())}")
    design.rebuild()

    print("[2/3] Importing SQDMetal PALACE eigenmode module...")
    from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
    from SQDMetal.Utilities.Materials import MaterialInterface
    import matplotlib
    matplotlib.use("Agg", force=True)

    print("[3/3] Running a VERY COARSE eigenmode simulation (smoke test only)...")
    t1 = time.time()
    user_defined_options = {
        "mesh_refinement": 1,
        "dielectric_material": "silicon",
        "starting_freq": 7.5e9,
        "number_of_freqs": 2,       # coarse: just enough to see something converge
        "solns_to_save": 2,
        "solver_order": 1,
        "solver_tol": 1.0e-5,       # loose tolerance, faster convergence
        "solver_maxits": 300,
        "fillet_resolution": 8,     # coarse fillet approximation
        "palace_dir": args.palace_bin,
        "num_cpus": args.num_cpus,
    }

    eigen_sim = PALACE_Eigenmode_Simulation(
        name="coarse_smoke_test",
        metal_design=design,
        sim_parent_directory=args.outdir,
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

    # Taper sizes informed by the 2026-08-25 AMR run: the Indicator field
    # showed AMR concentrating refinement around the readout resonator
    # (top-indicator points clustered along its trace edges, spanning a
    # wide area) and barely touching the qubit pad/SQUID arms (~1-4% of
    # total indicator mass there). So the readout's fine-mesh zone is
    # widened (taper_dist_min up) and the qubit's is narrowed (taper_dist_min
    # down) to match where refinement actually matters.
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

    eigen_sim.prepare_simulation()
    eigen_sim.run()
    print(f"      Palace run OK ({time.time()-t1:.1f}s)")

    port_epr = eigen_sim.retrieve_mode_port_EPR()
    mesh_sizes = eigen_sim.retrieve_simulation_sizes()
    print(f"      Mesh: {mesh_sizes}")
    print(f"      Eigenfrequencies (GHz): {[f/1e9 for f in port_epr['eigenfrequencies']]}")
    print(f"      Loaded Q: {list(port_epr['loaded_Q'])}")
    print(f"      mat_mode_port shape: "
          f"{None if port_epr['mat_mode_port'] is None else port_epr['mat_mode_port'].shape}")
    if port_epr['mat_mode_port'] is not None:
        print(f"      mat_mode_port:\n{port_epr['mat_mode_port']}")

    print(f"\n=== SMOKE TEST COMPLETE in {time.time()-t0:.1f}s ===")


if __name__ == "__main__":
    main()
