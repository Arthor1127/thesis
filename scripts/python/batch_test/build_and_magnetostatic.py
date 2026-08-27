"""
Magnetostatic (inductance-matrix) simulation, mirroring the structure of
build_and_eigenmode.py / build_and_driven.py.

Unlike the capacitance route, ONE Palace run produces the full NxN matrix:
PALACE_Inductance_Simulation collects every registered current source into
the config's "SurfaceCurrent" array, and Palace solves once per source
internally (unit current in source i, all others open).

Current sources are U-clips: a coplanar metallic piece that leaves the chip
edge and wraps around to join the two ground planes, closing the current
path. Without a closed path a magnetostatic solve is meaningless, so every
driven launcher/route MUST get a clip.

Output: retrieve_data() returns flux per unit current (Wb/A == H). Written
to inductance.json in HENRIES - note Chip.set_inductance_matrix() defaults
bare arrays to nH, so pass units="H" explicitly when ingesting.

Usage:
    python build_and_magnetostatic.py --ground-pos none \
        --outdir ~/ind_v1/geom_none --palace-bin $PALACE_BIN --num-cpus 16
"""
import argparse
import json
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import build_design, TOTAL_LENGTH_MM


def parse_ground_pos(s):
    if s.lower() == "none":
        return None
    return float(s)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ground-pos", type=parse_ground_pos, required=True,
                   help="mm, or the literal string 'none' for the ungrounded baseline")
    p.add_argument("--outdir", required=True)
    p.add_argument("--palace-bin", required=True)
    p.add_argument("--num-cpus", type=int, default=16)
    p.add_argument("--total-length-mm", type=float, default=None)
    p.add_argument("--launchers", nargs="+", default=["LP1", "LP2"],
                   help="LaunchpadWirebond components to drive. Each gets a U-clip "
                        "and becomes one row/column of the inductance matrix.")
    p.add_argument("--clip-side-um", type=float, default=20.0)
    p.add_argument("--clip-back-um", type=float, default=20.0)
    p.add_argument("--clip-gap-um", type=float, default=20.0)
    p.add_argument("--solver-order", type=int, default=2)
    p.add_argument("--bg-min-size-um", type=float, default=15.0)
    p.add_argument("--bg-max-size-um", type=float, default=150.0)
    args = p.parse_args()

    ground_pos_label = "none" if args.ground_pos is None else f"{args.ground_pos}mm"
    total_length_mm = args.total_length_mm if args.total_length_mm is not None else TOTAL_LENGTH_MM

    t0 = time.time()
    print(f"[1/5] Building design (ground_pos={ground_pos_label}, "
          f"total_length={total_length_mm}mm)...")
    design = build_design(ground_pos_mm=args.ground_pos, total_length_mm=total_length_mm)
    design.rebuild()
    print(f"      OK ({time.time()-t0:.1f}s). Components: {list(design.components.keys())}")

    print("[2/5] Importing SQDMetal PALACE inductance module...")
    from SQDMetal.PALACE.Inductance_Simulation import PALACE_Inductance_Simulation
    print("      OK")

    import matplotlib
    matplotlib.use("Agg", force=True)   # must come AFTER the SQDMetal import

    print("[3/5] Setting up magnetostatic simulation...")
    t1 = time.time()
    user_options = {
        # --- PALACE_Inductance_Simulation.default_user_options ---
        "solns_to_save": -1,            # -1 => one field solution per current source
        "solver_order": args.solver_order,
        "solver_tol": 1.0e-8,
        "solver_maxits": 100,
        # --- PALACE_Model_Base.default_user_options_parent ---
        "dielectric_material": "silicon",
        "fillet_resolution": 12,
        "mesh_min": args.bg_min_size_um * 1e-6,
        "mesh_max": args.bg_max_size_um * 1e-6,
        "palace_dir": args.palace_bin,
        "num_cpus": args.num_cpus,
    }

    try:
        # NOTE: signature is (name, sim_parent_directory, mode, meshing, ...).
        # There is NO metal_design positional - it goes through **kwargs.
        ind_sim = PALACE_Inductance_Simulation(
            name=f"ind_gpos_{ground_pos_label}",
            sim_parent_directory=args.outdir,
            mode="simPC",
            meshing="GMSH",
            user_options=user_options,
            create_files=True,
            metal_design=design,
        )
        ind_sim.add_metallic(1)
        ind_sim.add_ground_plane()

        for lp in args.launchers:
            ind_sim.create_current_source_with_Uclip_on_Launcher(
                lp,
                thickness_side=args.clip_side_um * 1e-6,
                thickness_back=args.clip_back_um * 1e-6,
                separation_gap=args.clip_gap_um * 1e-6,
            )
            print(f"      current source + U-clip on {lp}")

        ind_sim.fine_mesh_components(
            args.launchers + ["TL", "resonator1"],
            min_size=args.bg_min_size_um * 1e-6, max_size=args.bg_max_size_um * 1e-6,
            taper_dist_min=10e-6, metals_only=False,
        )

        print("[4/5] Meshing + running Palace...")
        ind_sim.prepare_simulation()
        ind_sim.run()
        print(f"      OK ({time.time()-t1:.1f}s)")

        print("[5/5] Retrieving inductance data...")
        data = ind_sim.retrieve_data()
        print(data)

        out = {
            "ground_pos_mm": args.ground_pos,
            "total_length_mm": total_length_mm,
            "launchers": args.launchers,
            "units": "H",
            "note": "flux per unit current (Wb/A). Chip.set_inductance_matrix "
                    "defaults to nH - pass units='H'.",
            "data": data if isinstance(data, (dict, list)) else str(data),
        }
        path = os.path.join(args.outdir, "inductance.json")
        with open(path, "w") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"      wrote {path}")

    except Exception:
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
