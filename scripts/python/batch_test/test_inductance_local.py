"""
Minimal LOCAL smoke test for PALACE_Inductance_Simulation.

Deliberately tiny: a 3x2 mm chip with a straight CPW between two edge
launchers. That is the smallest geometry that still exercises the full
path - U-clip current sources, two SurfaceCurrent entries, a 2x2
inductance matrix - without needing cluster memory.

The point is to answer three questions before spending cluster time:
  1. Does the class import and construct with this SQDMetal version?
  2. Do the U-clips render and produce valid current sources?
  3. What EXACTLY does retrieve_data() return?

A straight CPW also has a rough analytic sanity check: the self
inductance should land near L_l * length, with L_l ~ 0.4 uH/m for a
typical 10um/6um CPW on silicon. For a 2mm line that is ~0.8 nH. If you
get microhenries or picohenries, something is wrong with the geometry or
the units.

Usage:
    python test_inductance_local.py --palace-bin /path/to/palace
    python test_inductance_local.py --palace-bin $(which palace) --mesh-only
"""
import argparse
import json
import os
import sys
import time
import traceback


def build_test_design(chip_x_mm=3.0, chip_y_mm=2.0, launcher_offset_mm=1.35):
    """Two edge launchers joined by a straight CPW. Nothing meandered,
    nothing filleted into a corner - as boring as possible on purpose."""
    from qiskit_metal import designs, Dict
    from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
    from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder

    design = designs.DesignPlanar()
    design.overwrite_enabled = True
    design.chips.main.size.size_x = f"{chip_x_mm}mm"
    design.chips.main.size.size_y = f"{chip_y_mm}mm"

    LaunchpadWirebond(design, "LP1", options=Dict(
        pos_x=f"-{launcher_offset_mm}mm", pos_y="0mm", orientation="0",
        lead_length="30um", trace_width="10um", trace_gap="6um"))
    LaunchpadWirebond(design, "LP2", options=Dict(
        pos_x=f"{launcher_offset_mm}mm", pos_y="0mm", orientation="180",
        lead_length="30um", trace_width="10um", trace_gap="6um"))

    RoutePathfinder(design, "TL", options=Dict(
        trace_width="10um", trace_gap="6um",
        hfss_wire_bonds=False,
        pin_inputs=Dict(
            start_pin=Dict(component="LP1", pin="tie"),
            end_pin=Dict(component="LP2", pin="tie"))))

    design.rebuild()
    return design


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--palace-bin", default=os.environ.get("PALACE_BIN", "palace"),
                   help="Path to the palace executable. Defaults to $PALACE_BIN, then 'palace'.")
    p.add_argument("--palace-mode", default="local", choices=["local", "wsl", "docker"])
    p.add_argument("--outdir", default=os.path.expanduser("~/palace_ind_test"),
                   help="Parent directory for simulation files.")
    p.add_argument("--num-cpus", type=int, default=4)
    p.add_argument("--mesh-min-um", type=float, default=20.0)
    p.add_argument("--mesh-max-um", type=float, default=200.0)
    p.add_argument("--solver-order", type=int, default=1,
                   help="1 for a fast smoke test. Bump to 2 for anything you'd quote.")
    p.add_argument("--mesh-only", action="store_true",
                   help="Build + mesh + write the config, then stop. Verifies "
                        "geometry and the generated JSON without a solve.")
    p.add_argument("--view-gmsh", action="store_true",
                   help="Open the GMSH GUI to eyeball the U-clips. Needs a display.")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    if not args.outdir.endswith("/"):
        args.outdir += "/"

    t0 = time.time()
    print("[1/5] Building minimal test design...")
    design = build_test_design()
    print(f"      OK ({time.time()-t0:.1f}s). Components: {list(design.components.keys())}")

    print("[2/5] Importing PALACE_Inductance_Simulation...")
    from SQDMetal.PALACE.Inductance_Simulation import PALACE_Inductance_Simulation
    import matplotlib
    matplotlib.use("Agg", force=True)   # AFTER the SQDMetal import, not before
    print("      OK")

    user_options = {
        "solns_to_save": -1,
        "solver_order": args.solver_order,
        "solver_tol": 1.0e-8,
        "solver_maxits": 100,
        "dielectric_material": "silicon",
        "fillet_resolution": 12,
        "mesh_min": args.mesh_min_um * 1e-6,
        "mesh_max": args.mesh_max_um * 1e-6,
        "palace_dir": args.palace_bin,
        "palace_mode": args.palace_mode,
        "num_cpus": args.num_cpus,
    }

    try:
        print("[3/5] Constructing simulation object...")
        # Signature: (name, sim_parent_directory, mode, meshing, ...).
        # metal_design has NO positional slot - it rides in **kwargs.
        sim = PALACE_Inductance_Simulation(
            name="ind_smoketest",
            sim_parent_directory=args.outdir,
            mode="simPC",
            meshing="GMSH",
            user_options=user_options,
            view_design_gmsh_gui=args.view_gmsh,
            create_files=True,
            metal_design=design,
        )
        sim.add_metallic(1)
        sim.add_ground_plane()

        for lp in ["LP1", "LP2"]:
            sim.create_current_source_with_Uclip_on_Launcher(
                lp, thickness_side=20e-6, thickness_back=20e-6, separation_gap=20e-6)
            print(f"      current source + U-clip on {lp}")
        print(f"      registered {len(sim._ports)} current source(s) "
              f"-> expect a {len(sim._ports)}x{len(sim._ports)} matrix")
                # NOT optional: GMSH_Mesh_Builder.build_mesh() leaves thresh_field_id
        # unbound if no fine-mesh elements were registered.
        sim.fine_mesh_components(
            ["TL", "LP1", "LP2"],
            min_size=args.mesh_min_um * 1e-6,
            max_size=args.mesh_max_um * 1e-6,
            taper_dist_min=10e-6, metals_only=False,
        )
        
        print("[4/5] Meshing...")
        t1 = time.time()
        sim.prepare_simulation()
        print(f"      OK ({time.time()-t1:.1f}s)")

        cfg = os.path.join(args.outdir, "ind_smoketest", "ind_smoketest.json")
        if os.path.exists(cfg):
            with open(cfg) as f:
                conf = json.load(f)
            sc = conf.get("Boundaries", {}).get("SurfaceCurrent", [])
            print(f"      config: Problem.Type={conf.get('Problem', {}).get('Type')}, "
                  f"{len(sc)} SurfaceCurrent boundaries")
            if conf.get("Problem", {}).get("Type") != "Magnetostatic":
                print("      !!! Problem.Type is not Magnetostatic - stop and investigate")

        if args.mesh_only:
            print("[5/5] --mesh-only given, stopping before the solve.")
            print(f"      Inspect: {cfg}")
            return

        print("[5/5] Running Palace...")
        t2 = time.time()
        sim.run()
        print(f"      solve done ({time.time()-t2:.1f}s)")

        data = sim.retrieve_data()
        print("\n" + "=" * 60)
        print("retrieve_data() returned:")
        print(f"  type: {type(data)}")
        print(data)
        print("=" * 60)
        print("Units are Wb/A == henries. A ~2mm CPW should be order 1 nH.")

        with open(os.path.join(args.outdir, "smoketest_result.json"), "w") as f:
            json.dump({"raw": data if isinstance(data, (dict, list)) else str(data)},
                      f, indent=2, default=str)

    except Exception:
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
