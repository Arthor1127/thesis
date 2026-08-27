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
    p.add_argument("--ground-pos", type=parse_ground_pos, required=True)
    p.add_argument("--freq-start-ghz", type=float, default=14.0,
                    help="Sweep start (GHz). Was a hardcoded module constant.")
    p.add_argument("--freq-end-ghz", type=float, default=16.0,
                    help="Sweep end (GHz).")
    p.add_argument("--n-points", type=int, default=2001,
                    help="Number of frequency points REPORTED. With "
                         "--adaptive-tol > 0 this is NOT the number of "
                         "full-order solves - Palace solves at a few "
                         "auto-chosen samples (<= --adaptive-max-samples) and "
                         "evaluates the rest from the reduced-order model, so "
                         "a dense grid here is nearly free. That matters here: "
                         "the eigenmode Q of ~4e5 implies an internal "
                         "linewidth of ~37 kHz, which no viable uniform grid "
                         "over a 2 GHz span could resolve (it would need "
                         "~500k points). The loaded Q with ports attached will "
                         "be lower - possibly much lower, since coupling to "
                         "the 50 ohm feedline usually dominates - but the ROM "
                         "is what makes a sharp response representable at all.")
    p.add_argument("--outdir", required=True)
    p.add_argument("--palace-bin", required=True)
    p.add_argument("--num-cpus", type=int, default=16)
    p.add_argument("--total-length-mm", type=float, default=None)
    p.add_argument("--bg-min-size-um", type=float, default=15.0)
    p.add_argument("--bg-max-size-um", type=float, default=150.0)
    p.add_argument("--resonator-bg-min-size-um", type=float, default=None)
    p.add_argument("--resonator-bg-max-size-um", type=float, default=None)
    p.add_argument("--chip-size-x-mm", type=float, default=None)
    p.add_argument("--chip-size-y-mm", type=float, default=None)
    p.add_argument("--chip-center-x-mm", type=float, default=None)
    p.add_argument("--chip-center-y-mm", type=float, default=None)
    p.add_argument("--strap-margin-um", type=float, default=3.0)
    p.add_argument("--strap-taper-dist-min-um", type=float, default=2.0,
                    help="Within this distance of the strap box, the full "
                         "--resonator-min-size-um applies.")
    p.add_argument("--strap-taper-dist-max-um", type=float, default=40.0,
                    help="Distance over which the strap box's fine mesh "
                         "relaxes back to the resonator background. Was "
                         "hardcoded at 6um, which compresses a ~75-150x size "
                         "ratio into a 4um band - steep enough that it "
                         "destabilised AMR on the eigenmode side until "
                         "widened. Defaulted wider here for the same reason.")
    p.add_argument("--solver-order", type=int, default=1)
    p.add_argument("--adaptive-tol", type=float, default=1e-3)
    p.add_argument("--adaptive-max-samples", type=int, default=30)
    p.add_argument("--solns-to-save", type=int, default=10)
    p.add_argument("--solver-tol", type=float, default=1.0e-6)
    p.add_argument("--resonator-min-size-um", type=float, default=2.0)
    args = p.parse_args()

    resonator_bg_min_um = (args.resonator_bg_min_size_um
                            if args.resonator_bg_min_size_um is not None
                            else args.bg_min_size_um)
    resonator_bg_max_um = (args.resonator_bg_max_size_um
                            if args.resonator_bg_max_size_um is not None
                            else args.bg_max_size_um)

    ground_pos_label = "none" if args.ground_pos is None else f"{args.ground_pos}mm"
    total_length_mm = args.total_length_mm if args.total_length_mm is not None else TOTAL_LENGTH_MM
    t0 = time.time()

    freq_start_hz = args.freq_start_ghz * 1e9
    freq_end_hz = args.freq_end_ghz * 1e9
    freq_step = (freq_end_hz - freq_start_hz) / max(args.n_points - 1, 1)
    print(f"[1/5] Continuous sweep: {args.freq_start_ghz:.4f} - {args.freq_end_ghz:.4f} GHz "
          f"({args.n_points} points, step={freq_step/1e6:.4f} MHz)")
    print(f"      Mesh background: ports(TL/LP1/LP2)={args.bg_min_size_um}/"
          f"{args.bg_max_size_um}um, resonator1={resonator_bg_min_um}/"
          f"{resonator_bg_max_um}um, strap-local={args.resonator_min_size_um}um")
    print(f"      Solver: order={args.solver_order}, tol={args.solver_tol:.1e}, "
          f"solns_to_save={args.solns_to_save}")

    print(f"[2/5] Building design (ground_pos={ground_pos_label}, "
          f"total_length={total_length_mm}mm)...")
    design = build_design(ground_pos_mm=args.ground_pos,
                            total_length_mm=total_length_mm,
                            chip_size_x_mm=args.chip_size_x_mm,
                            chip_size_y_mm=args.chip_size_y_mm,
                            chip_center_x_mm=args.chip_center_x_mm,
                            chip_center_y_mm=args.chip_center_y_mm)
    print(f"      OK ({time.time()-t0:.1f}s)")

    # Ensure qgeometry tables (needed below to locate the strap, if any)
    # are populated before we read them.
    design.rebuild()

    print("[3/5] Importing SQDMetal PALACE driven-simulation module...")
    from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation
    print("      OK")

    import matplotlib
    matplotlib.use("Agg", force=True)

    print("[4/5] Setting up + running driven simulation (production settings)...")
    t1 = time.time()
    driven_options = {
        "dielectric_material": "silicon",
        "solns_to_save": args.solns_to_save,
        "solver_order": args.solver_order,
        "solver_tol": args.solver_tol,
        "solver_maxits": 500,
        "solver_initial_guess": True,  # warm-start each freq step from the previous one
        "fillet_resolution": 12,
        "palace_dir": args.palace_bin,
        "num_cpus": args.num_cpus,
    }

    try:
        driven_sim = PALACE_Driven_Simulation(
            name=f"driven_gpos_{ground_pos_label}",
            metal_design=design,
            sim_parent_directory=args.outdir,
            mode="simPC",
            meshing="GMSH",
            user_options=driven_options,
            create_files=True,
        )
        driven_sim.add_metallic(1)
        driven_sim.add_ground_plane()
        driven_sim.create_port_CPW_on_Launcher("LP1", 20e-6)
        driven_sim.create_port_CPW_on_Launcher("LP2", 20e-6)
        driven_sim.set_port_excitation(1)

        driven_sim.fine_mesh_components(
            ["TL", "LP1", "LP2"],
            min_size=args.bg_min_size_um * 1e-6, max_size=args.bg_max_size_um * 1e-6,
            taper_dist_min=10e-6, metals_only=False,
        )
        driven_sim.fine_mesh_components(
            ["resonator1"],
            min_size=resonator_bg_min_um * 1e-6, max_size=resonator_bg_max_um * 1e-6,
            taper_dist_min=10e-6, metals_only=False,
        )

        if args.ground_pos is not None:
            poly_table = design.qgeometry.tables["poly"]
            strap_rows = poly_table[poly_table["name"].str.startswith("ground_strap_")]
            if len(strap_rows) == 0:
                print("WARNING: --ground-pos was given but no ground_strap_* "
                      "geometry was found in the design - the strap-local fine "
                      "mesh will not be applied.")
            else:
                strap_geom = strap_rows.iloc[0]["geometry"]
                minx, miny, maxx, maxy = strap_geom.bounds
                safe_margin_mm = args.strap_margin_um / 1000

                x1, y1 = (minx - safe_margin_mm) * 1e-3, (miny - safe_margin_mm) * 1e-3
                x2, y2 = (maxx + safe_margin_mm) * 1e-3, (maxy + safe_margin_mm) * 1e-3
                print(f"      Strap bounding box (mm): x=[{minx:.4f},{maxx:.4f}] "
                      f"y=[{miny:.4f},{maxy:.4f}], margin={safe_margin_mm*1000:.1f}um - "
                      f"applying {args.resonator_min_size_um}um fine mesh there only.")

                driven_sim.fine_mesh_in_rectangle(
                    x1, y1, x2, y2,
                    min_size=args.resonator_min_size_um * 1e-6, max_size=resonator_bg_max_um * 1e-6,
                    taper_dist_min=args.strap_taper_dist_min_um * 1e-6,
                    taper_dist_max=args.strap_taper_dist_max_um * 1e-6,
                )

        driven_sim.set_freq_values(
            freq_start=freq_start_hz,
            freq_end=freq_end_hz,
            freq_step=freq_step,
        )

        driven_sim.prepare_simulation()

        if args.adaptive_tol > 0:
            cfg_path = driven_sim._sim_config
            with open(cfg_path, "r") as f:
                cfg = json.load(f)
            cfg["Solver"]["Driven"]["AdaptiveTol"] = args.adaptive_tol
            cfg["Solver"]["Driven"]["AdaptiveMaxSamples"] = args.adaptive_max_samples
            with open(cfg_path, "w") as f:
                json.dump(cfg, f, indent=2)
            print(f"      Adaptive fast frequency sweep ENABLED: "
                  f"AdaptiveTol={args.adaptive_tol:.1e}, "
                  f"AdaptiveMaxSamples={args.adaptive_max_samples} "
                  f"(vs {args.n_points} full-order solves if disabled)")
        else:
            print(f"      Adaptive sweep disabled - solving full-order model "
                  f"at all {args.n_points} frequencies (slow; pass "
                  f"--adaptive-tol 1e-3 to enable PROM)")

        driven_sim.run()
    except Exception:
        print("\n!!! FAILED during driven simulation !!!")
        traceback.print_exc()
        sys.exit(1)

    print(f"      OK ({time.time()-t1:.1f}s)")

    print("[5/5] Retrieving S-parameter data...")
    s_data = driven_sim.retrieve_data()
    # SQDMetal returns a dict keyed 'freqs' (plural). An earlier version of
    # this script checked 'freq' (singular), which reported a genuinely
    # successful Palace run as a failure.
    if s_data is None or "freqs" not in s_data:
        print("\n!!! FAILED: run() completed without raising, but produced no "
              "usable output (retrieve_data() returned None or incomplete data).")
        print("This usually means the MPI job itself was killed mid-run "
              "(check for OOM/signal 9 above) without Palace/SQDMetal "
              "propagating that as a Python exception.")
        sys.exit(1)

    # Persist the result. Without this the sweep produces nothing an
    # aggregator can read - Palace's own port-S.csv survives inside the
    # simulation directory, but each array task needs its S21 in a single
    # predictable place keyed by strap position.
    import numpy as np
    freqs = np.asarray(s_data["freqs"])
    s11 = np.asarray(s_data.get("S11"))
    s21 = np.asarray(s_data.get("S21"))

    npz_path = f"{args.outdir}/s_sweep.npz"
    np.savez(npz_path,
             ground_pos_mm=(-1.0 if args.ground_pos is None else args.ground_pos),
             ground_pos_is_baseline=(args.ground_pos is None),
             total_length_mm=total_length_mm,
             freq_hz=freqs, S11=s11, S21=s21)

    # Locate the |S21| minimum - for a notch-type response through a
    # feedline this is where the resonance sits. Reported as a quick sanity
    # check that the swept band actually brackets a mode; if the minimum
    # lands on the first or last point, the band is probably wrong.
    s21_db = 20 * np.log10(np.abs(s21)) if s21.size else np.array([])
    summary = {
        "status": "success",
        "ground_pos_mm": args.ground_pos,
        "total_length_mm": total_length_mm,
        "freq_start_ghz": args.freq_start_ghz,
        "freq_end_ghz": args.freq_end_ghz,
        "n_points": int(freqs.size),
        "solver_order": args.solver_order,
        "adaptive_tol": args.adaptive_tol,
        "wall_time_s": time.time() - t0,
    }
    if s21_db.size:
        i = int(np.argmin(s21_db))
        summary["s21_min_db"] = float(s21_db[i])
        summary["s21_min_freq_ghz"] = float(freqs[i] / 1e9)
        edge = (i == 0 or i == freqs.size - 1)
        summary["s21_min_at_band_edge"] = edge
        print(f"      |S21| minimum: {s21_db[i]:.2f} dB at "
              f"{freqs[i]/1e9:.6f} GHz")
        if edge:
            print("      WARNING: the minimum sits at a band edge, so the "
                  "resonance is probably outside "
                  f"{args.freq_start_ghz}-{args.freq_end_ghz} GHz. Widen the "
                  "band or re-check the expected mode frequency.")

    json.dump(summary, open(f"{args.outdir}/sweep_summary.json", "w"), indent=2)
    print(f"      Saved {npz_path} ({freqs.size} points)")
    print(f"\n=== DRIVEN SWEEP COMPLETE in {summary['wall_time_s']:.1f}s ===")


if __name__ == "__main__":
    main()