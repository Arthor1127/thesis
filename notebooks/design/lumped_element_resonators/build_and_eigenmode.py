import argparse
import json
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import (build_design, CAP_FINGER_COUNT, CAP_GAP_UM,
                           RES_TOTAL_LENGTH_MM)


def _json_safe(obj):
    """Recursively convert numpy arrays/scalars into plain Python types so
    json.dump doesn't choke on them.

    numpy complex scalars have .tolist(), but it converts them to a plain
    Python `complex` - which is STILL not JSON serializable, so that result
    must be re-checked (see build_and_eigenmode.py in squid_resonator_epr,
    which hit this exact issue on interface_epr data)."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, complex):
        return {"real": obj.real, "imag": obj.imag}
    if hasattr(obj, "tolist"):
        return _json_safe(obj.tolist())
    return obj


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--finger-count", type=int, default=None,
                    help="Number of capacitor fingers. Defaults to build_design.CAP_FINGER_COUNT.")
    p.add_argument("--cap-gap-um", type=float, default=None,
                    help="Capacitor finger gap (um). Defaults to build_design.CAP_GAP_UM.")
    p.add_argument("--res-length-mm", type=float, default=None,
                    help="Meander inductor total length (mm). Defaults to build_design.RES_TOTAL_LENGTH_MM.")
    p.add_argument("--outdir", required=True)
    p.add_argument("--palace-bin", required=True)
    p.add_argument("--num-cpus", type=int, default=16)
    p.add_argument("--number-of-freqs", type=int, default=3,
                    help="Number of eigenmodes to solve. A bare LC loop only "
                        "has one mode of interest (the fundamental) plus its "
                        "harmonics, so this can be much smaller than for a "
                        "qubit+readout design.")
    p.add_argument("--starting-freq", type=float, default=1e9)

    # SQDMetal's mesh sizing is a single Gmsh "Threshold" field per
    # fine_mesh_components(...) call: SizeMin within taper_dist_min of the
    # cap+inductor geometry, ramping linearly up to SizeMax by taper_dist_max,
    # then FLAT AT SizeMax everywhere beyond that - all the way to the chip
    # edges (see GMSH_Mesh_Builder.build_mesh(): the fields list is built
    # only from explicit fine_mesh_components calls, there is no separate
    # implicit "background" field to combine against). Since cap+inductor is
    # the entire metal footprint of this design, --max-size-um IS the mesh
    # size for the rest of the chip too, not just "near the resonator" - an
    # earlier version of this script exposed separate --bg-*/--resonator-*
    # flags as if they were independent, but the --bg-* ones never reached
    # Gmsh and changing them had no effect. There is only one real lever.
    p.add_argument("--min-size-um", type=float, default=6.0,
                    help="Mesh element size (um) right at the cap+inductor metal. "
                        "6um resolves cap_gap=6um / cpw_gap by default.")
    p.add_argument("--max-size-um", type=float, default=100.0,
                    help="Mesh element size (um) once --taper-dist-max-um from the metal is "
                        "reached - this is also the effective size for the rest of the chip "
                        "(ground plane, substrate, vacuum), since no other field exists.")
    p.add_argument("--taper-dist-min-um", type=float, default=40.0,
                    help="Distance (um) over which mesh elements stay at --min-size-um "
                        "before beginning to relax outward.")
    p.add_argument("--taper-dist-max-um", type=float, default=200.0,
                    help="Distance (um) at which mesh elements reach --max-size-um and "
                        "flatten out for the remainder of the chip.")

    p.add_argument("--amr-max-its", type=int, default=0)
    p.add_argument("--amr-tol", type=float, default=1e-2)
    p.add_argument("--chip-size-x-mm", type=float, default=None)
    p.add_argument("--chip-size-y-mm", type=float, default=None)
    p.add_argument("--chip-center-x-mm", type=float, default=None)
    p.add_argument("--chip-center-y-mm", type=float, default=None)
    args = p.parse_args()

    finger_count = args.finger_count if args.finger_count is not None else CAP_FINGER_COUNT
    cap_gap_um = args.cap_gap_um if args.cap_gap_um is not None else CAP_GAP_UM
    res_length_mm = args.res_length_mm if args.res_length_mm is not None else RES_TOTAL_LENGTH_MM

    t0 = time.time()
    print(f"[1/4] Building design (finger_count={finger_count}, cap_gap={cap_gap_um}um, "
          f"res_length={res_length_mm}mm)...")
    design = build_design(finger_count=finger_count, cap_gap_um=cap_gap_um,
                            res_length_mm=res_length_mm,
                            chip_size_x_mm=args.chip_size_x_mm,
                            chip_size_y_mm=args.chip_size_y_mm,
                            chip_center_x_mm=args.chip_center_x_mm,
                            chip_center_y_mm=args.chip_center_y_mm)
    print(f"      OK ({time.time()-t0:.1f}s). Components: {list(design.components.keys())}")

    design.rebuild()

    print("[2/4] Importing SQDMetal PALACE eigenmode module...")
    from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
    from SQDMetal.Utilities.Materials import MaterialInterface
    print("      OK")

    import matplotlib
    matplotlib.use("Agg", force=True)

    print("[3/4] Setting up + running eigenmode simulation...")
    t1 = time.time()
    user_defined_options = {
        "mesh_refinement": 0,  # UniformLevels - NOT AMR; see --amr-max-its
        "dielectric_material": "silicon",
        "starting_freq": args.starting_freq,
        "number_of_freqs": args.number_of_freqs,
        "solns_to_save": args.number_of_freqs,
        "solver_order": 1,
        "solver_tol": 1.0e-5,
        "solver_maxits": 300,
        "fillet_resolution": 12,
        "palace_dir": args.palace_bin,
        "num_cpus": args.num_cpus,
    }

    run_label = f"fingers_{finger_count}_capgap_{cap_gap_um}um_reslen_{res_length_mm}mm"

    def _fail(status, extra=None):
        payload = {"status": status, "finger_count": finger_count,
                   "cap_gap_um": cap_gap_um, "res_length_mm": res_length_mm}
        if extra:
            payload.update(extra)
        json.dump(payload, open(f"{args.outdir}/eigenmode_result.json", "w"), indent=2)

    try:
        eigen_sim = PALACE_Eigenmode_Simulation(
            name=f"eig_{run_label}",
            metal_design=design,
            sim_parent_directory=args.outdir,
            mode="simPC",
            meshing="GMSH",
            user_options=user_defined_options,
            create_files=True,
        )
        eigen_sim.add_metallic(1)
        eigen_sim.add_ground_plane()

        # No junctions/ports here - this is a bare LC loop (cap + meander),
        # not a qubit, so there is nothing to feed create_port_JosephsonJunction.
        # Eigenfrequencies come straight from eig.csv via retrieve_data().

        eigen_sim.fine_mesh_components(
            ["cap", "inductor"],
            min_size=args.min_size_um * 1e-6,
            max_size=args.max_size_um * 1e-6,
            taper_dist_min=args.taper_dist_min_um * 1e-6,
            taper_dist_max=args.taper_dist_max_um * 1e-6,
            metals_only=False,
        )

        eigen_sim.setup_EPR_interfaces(
            metal_air=MaterialInterface("Aluminium-Vacuum"),
            substrate_air=MaterialInterface("Silicon-Vacuum"),
            substrate_metal=MaterialInterface("Silicon-Aluminium"),
        )
        if args.amr_max_its > 0:
            eigen_sim._mesh_refinement = {
                "UniformLevels": 0,
                "MaxIts": args.amr_max_its,
                "Tol": args.amr_tol,
            }
            print(f"      Solution-based AMR ENABLED: MaxIts={args.amr_max_its}, "
                  f"Tol={args.amr_tol:.1e}")

        eigen_sim.prepare_simulation()
        eigen_sim.run()
    except Exception:
        print("\n!!! FAILED during eigenmode simulation !!!")
        traceback.print_exc()
        _fail("failed")
        sys.exit(1)

    print(f"      OK ({time.time()-t1:.1f}s)")

    print("[4/4] Retrieving eigenmode data and writing eigenmode_result.json...")
    eig_csv_path = eigen_sim._output_data_dir + "/eig.csv"
    if not os.path.exists(eig_csv_path):
        print(f"\n!!! FAILED: run() completed without raising, but {eig_csv_path} "
              "was never written.")
        print("This usually means the MPI job itself was killed mid-run "
              "(check for OOM/signal 9 in the log above) without Palace/SQDMetal "
              "propagating that as a Python exception.")
        _fail("failed", {"reason": "missing eig.csv, likely killed MPI job"})
        sys.exit(1)

    import pandas as pd
    eig_data = pd.read_csv(eig_csv_path)
    eig_data.columns = [c.strip() for c in eig_data.columns]

    try:
        interface_epr = eigen_sim.retrieve_interface_EPR_data()
    except Exception:
        print("      WARNING: retrieve_interface_EPR_data() failed - continuing without it.")
        traceback.print_exc()
        interface_epr = None

    mesh_sizes = eigen_sim.retrieve_simulation_sizes()

    num_modes = len(eig_data)
    if num_modes == 0:
        print(f"\n!!! FAILED: zero converged modes. starting_freq="
              f"{args.starting_freq/1e9:.3f}GHz, number_of_freqs={args.number_of_freqs}.")
        _fail("no_modes_converged", {"starting_freq_hz": args.starting_freq,
                                      "number_of_freqs": args.number_of_freqs})
        sys.exit(1)

    modes = []
    for _, row in eig_data.iterrows():
        modes.append({
            "m": int(row["m"]),
            "f_ghz": float(row["Re{f} (GHz)"]),
            "Q": float(row["Q"]),
        })

    print("      Modes found:")
    for mode in modes:
        print(f"        m={mode['m']}: f={mode['f_ghz']:.4f}GHz, Q={mode['Q']:.3e}")

    result = {
        "status": "success",
        "finger_count": finger_count,
        "cap_gap_um": cap_gap_um,
        "res_length_mm": res_length_mm,
        "modes": modes,
        "mesh": {
            "elements": mesh_sizes.get("MeshElements"),
            "min_size_um": args.min_size_um,
            "max_size_um": args.max_size_um,
            "taper_dist_min_um": args.taper_dist_min_um,
            "taper_dist_max_um": args.taper_dist_max_um,
        },
        "interface_epr": _json_safe(interface_epr),
        "wall_time_s": time.time() - t0,
    }
    out_path = f"{args.outdir}/eigenmode_result.json"
    json.dump(result, open(out_path, "w"), indent=2)

    print(f"\n=== EIGENMODE RUN COMPLETE in {result['wall_time_s']:.1f}s ===")
    print(f"Result written to {out_path}")


if __name__ == "__main__":
    main()
