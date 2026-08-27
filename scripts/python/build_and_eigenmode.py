"""
Stage A of the sweep: for one grounding-strap position (or the ungrounded
baseline, ground_pos_mm=None), build the design and run an eigenmode
simulation. Writes resonance.json - primarily useful as a validation/
reference dataset (the harmonic driven sweep in Stage B uses a FIXED
frequency grid around n*5GHz for every geometry point, not this per-point
resonance, so that mode-splitting/shifting is visible in the final
colormap relative to a common reference grid).

The ground_pos_mm=None run is the most important one to check first. It is
the ungrounded baseline, and its true f1 tells you whether TOTAL_LENGTH_MM
in build_design.py is correctly calibrated for a 5 GHz fundamental.

Usage:
    python build_and_eigenmode.py --ground-pos none \
        --outdir ~/sweep_v2/geom_none --palace-bin $PALACE_BIN --num-cpus 16
    python build_and_eigenmode.py --ground-pos 5.9 \
        --outdir ~/sweep_v2/geom_5.9 --palace-bin $PALACE_BIN --num-cpus 16
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
    p.add_argument("--mode-index", type=int, default=0)
    p.add_argument("--number-of-freqs", type=int, default=2)
    p.add_argument("--starting-freq", type=float, default=5e9,
                    help="Hz - SLEPc shift-invert search center; target fundamental is ~5 GHz")
    p.add_argument("--total-length-mm", type=float, default=None,
                    help="Override TOTAL_LENGTH_MM from build_design.py, for calibration sweeps. "
                         "If omitted, uses the module default.")
    p.add_argument("--bg-min-size-um", type=float, default=15.0,
                    help="Background mesh min_size (um) for resonator1, away from "
                         "the grounding strap. Notebook mesh-only testing showed "
                         "8um->126.6M elements, 15um->19.7M elements (84%% reduction) "
                         "with fillet_resolution making no measurable difference.")
    p.add_argument("--bg-max-size-um", type=float, default=150.0,
                    help="Background mesh max_size (um).")
    p.add_argument("--strap-margin-um", type=float, default=3.0,
                    help="Fixed padding (um) around the strap's bounding box. "
                         "Must stay above 0 - see build_and_driven.py for why a "
                         "clearance-derived version of this caused a crash.")
    args = p.parse_args()

    ground_pos_label = "none" if args.ground_pos is None else f"{args.ground_pos}mm"
    total_length_mm = args.total_length_mm if args.total_length_mm is not None else TOTAL_LENGTH_MM
    t0 = time.time()
    print(f"[1/4] Building design (ground_pos={ground_pos_label}, "
          f"total_length={total_length_mm}mm)...")
    design = build_design(ground_pos_mm=args.ground_pos, total_length_mm=total_length_mm)
    print(f"      OK ({time.time()-t0:.1f}s). Components: {list(design.components.keys())}")

    # Ensure qgeometry tables (needed below to locate the strap, if any)
    # are populated before we read them.
    design.rebuild()

    print("[2/4] Importing SQDMetal PALACE eigenmode module...")
    from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
    from SQDMetal.Utilities.Materials import MaterialInterface
    print("      OK")

    import matplotlib
    matplotlib.use("Agg", force=True)

    print("[3/4] Setting up + running eigenmode simulation (production settings)...")
    t1 = time.time()
    user_defined_options = {
        "mesh_refinement": 0,
        "dielectric_material": "silicon",
        "starting_freq": args.starting_freq,
        "number_of_freqs": args.number_of_freqs,
        "solns_to_save": args.number_of_freqs,
        "solver_order": 2,
        "solver_tol": 1.0e-8,
        "solver_maxits": 200,
        "fillet_resolution": 12,
        "palace_dir": args.palace_bin,
        "num_cpus": args.num_cpus,
    }

    try:
        eigen_sim = PALACE_Eigenmode_Simulation(
            name=f"eig_gpos_{ground_pos_label}",
            metal_design=design,
            sim_parent_directory=args.outdir,
            mode="simPC",
            meshing="GMSH",
            user_options=user_defined_options,
            create_files=True,
        )
        eigen_sim.add_metallic(1)
        eigen_sim.add_ground_plane()

        # Whole resonator at a coarse setting; a small box around the
        # actual strap (if any) gets the fine setting. This mirrors the
        # fix applied in build_and_driven.py - fine_mesh_components() is
        # component-scoped (the whole ~11.8mm resonator, not just the
        # ~24x17um strap), which was applying unnecessary fine resolution
        # across the entire trace, including on the ungrounded baseline
        # where there is no strap at all.
        eigen_sim.fine_mesh_components(
            ["resonator1"], min_size=args.bg_min_size_um * 1e-6, max_size=args.bg_max_size_um * 1e-6,
            taper_dist_min=10e-6, metals_only=False,
        )
        if args.ground_pos is not None:
            poly_table = design.qgeometry.tables["poly"]
            strap_rows = poly_table[poly_table["name"].str.startswith("ground_strap_")]
            if len(strap_rows) == 0:
                print("WARNING: --ground-pos was given but no ground_strap_* "
                      "geometry was found - the strap-local fine mesh will "
                      "not be applied.")
            else:
                strap_geom = strap_rows.iloc[0]["geometry"]
                minx, miny, maxx, maxy = strap_geom.bounds

                # Fixed margin (NOT clearance-derived). A previous version
                # measured distance from the strap to the resonator's own
                # trace and used a fraction of that as the margin - broken,
                # since the strap is physically attached to the trace at its
                # position (distance ~0 by design there). When that measured
                # distance came out as exactly 0, the margin collapsed to 0,
                # producing a degenerate fine-mesh region that crashed
                # Palace almost instantly. See build_and_driven.py for the
                # full writeup - same fix applied here.
                safe_margin_mm = args.strap_margin_um / 1000

                x1, y1 = (minx - safe_margin_mm) * 1e-3, (miny - safe_margin_mm) * 1e-3
                x2, y2 = (maxx + safe_margin_mm) * 1e-3, (maxy + safe_margin_mm) * 1e-3
                print(f"      Strap bounding box (mm): x=[{minx:.4f},{maxx:.4f}] "
                      f"y=[{miny:.4f},{maxy:.4f}], margin={safe_margin_mm*1000:.1f}um - "
                      f"applying 2um fine mesh there only.")
                eigen_sim.fine_mesh_in_rectangle(
                    x1, y1, x2, y2,
                    min_size=2e-6, max_size=args.bg_min_size_um * 1e-6,
                    taper_dist_min=2e-6, taper_dist_max=6e-6,
                )
        eigen_sim.setup_EPR_interfaces(
            metal_air=MaterialInterface("Aluminium-Vacuum"),
            substrate_air=MaterialInterface("Silicon-Vacuum"),
            substrate_metal=MaterialInterface("Silicon-Aluminium"),
        )
        eigen_sim.prepare_simulation()
        eigen_sim.run()
    except Exception:
        print("\n!!! FAILED during eigenmode simulation !!!")
        traceback.print_exc()
        json.dump({"status": "failed", "ground_pos_mm": args.ground_pos},
                   open(f"{args.outdir}/resonance.json", "w"), indent=2)
        sys.exit(1)

    print(f"      OK ({time.time()-t1:.1f}s)")

    print("[4/4] Parsing eig.csv output...")
    import os as _os
    eig_csv_path = eigen_sim._output_data_dir + "/eig.csv"
    if not _os.path.exists(eig_csv_path):
        print(f"\n!!! FAILED: run() completed without raising, but {eig_csv_path} "
              "was never written.")
        print("This usually means the MPI job itself was killed mid-run "
              "(check for OOM/signal 9 in the log above) without Palace/SQDMetal "
              "propagating that as a Python exception.")
        json.dump({"status": "failed", "ground_pos_mm": args.ground_pos,
                   "reason": "missing eig.csv, likely killed MPI job"},
                   open(f"{args.outdir}/resonance.json", "w"), indent=2)
        sys.exit(1)

    import pandas as pd
    eig_df = pd.read_csv(eig_csv_path)
    eig_df.columns = eig_df.columns.str.strip()
    f_re_col = [c for c in eig_df.columns if c.startswith("Re{f}")][0]
    q_col = [c for c in eig_df.columns if c.startswith("Q")][0]

    all_modes_ghz = eig_df[f_re_col].astype(float).tolist()
    resonance_ghz = float(eig_df.iloc[args.mode_index][f_re_col])
    q_factor = float(eig_df.iloc[args.mode_index][q_col])
    print(f"      All modes found (GHz): {all_modes_ghz}")
    print(f"      Mode {args.mode_index}: f = {resonance_ghz:.4f} GHz, Q = {q_factor:.3e}")

    if args.ground_pos is None:
        target = 5.0
        pct_off = 100 * (resonance_ghz - target) / target
        print(f"\n      *** BASELINE CHECK: measured f1={resonance_ghz:.4f} GHz vs "
              f"target 5.0000 GHz ({pct_off:+.2f}%) at total_length={total_length_mm}mm ***")
        if abs(pct_off) > 5:
            suggested = total_length_mm * resonance_ghz / target
            print(f"      WARNING: baseline is >5% off target. Suggested next length "
                  f"to try: --total-length-mm {suggested:.4f}")

    result = {
        "status": "success",
        "ground_pos_mm": args.ground_pos,
        "total_length_mm": total_length_mm,
        "mode_index": args.mode_index,
        "resonance_ghz": resonance_ghz,
        "q_factor": q_factor,
        "all_modes_ghz": all_modes_ghz,
        "wall_time_s": time.time() - t0,
    }
    out_path = f"{args.outdir}/resonance.json"
    json.dump(result, open(out_path, "w"), indent=2)

    print(f"\n=== STAGE A COMPLETE in {result['wall_time_s']:.1f}s ===")
    print(f"Result written to {out_path}")


if __name__ == "__main__":
    main()
