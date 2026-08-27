"""
Stage B of the sweep: for one grounding-strap position (or the ungrounded
baseline), run ONE continuous driven S-parameter sweep covering all four
harmonics simultaneously (a little below the 1st harmonic to a little
above the 4th), on a FIXED frequency grid shared by every geometry point.

Using the same grid for every position is what lets you see how S21 peaks
shift, split, or move as the grounding position changes, in the final
colormap (x=position, y=frequency, z=magnitude of S21).

--- On point density ---
The eigenmode simulation measured an internal (unloaded) Q of about 4e5.
At that Q, a resonance linewidth is f/Q, which at 5 GHz is only about
12.5 kHz - spanning a 16 GHz range with that resolution would need well
over a million points, which is not computationally feasible for a driven
FEM sweep. In practice, the driven simulation loaded Q (which includes
the 50 ohm port coupling) will likely be substantially lower and the
linewidth correspondingly broader, but that is not known in advance. The
default n_points below (321, about 50 MHz spacing across 16 GHz) is a
starting point, not a guarantee. If the resulting S21 plot shows no
visible peaks, the loaded linewidth is narrower than this spacing and
n_points needs to be increased, or a finer sweep run only near the actual
peak locations once they are known.

Usage (single position, for testing):
    python build_and_driven.py --ground-pos 5.906 \
        --outdir ~/sweep_v2/driven_5.906 --palace-bin $PALACE_BIN --num-cpus 16
    python build_and_driven.py --ground-pos none \
        --outdir ~/sweep_v2/driven_none --palace-bin $PALACE_BIN --num-cpus 16
"""
import argparse
import json
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import build_design, TOTAL_LENGTH_MM

FREQ_START_HZ = 4.5e9   # a little below the 1st harmonic (~5 GHz)
FREQ_END_HZ = 20.5e9    # a little above the 4th harmonic (~20 GHz)


def parse_ground_pos(s):
    if s.lower() == "none":
        return None
    return float(s)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ground-pos", type=parse_ground_pos, required=True)
    p.add_argument("--n-points", type=int, default=321,
                    help="points across the full 4.5-20.5 GHz range, ~50 MHz "
                         "spacing at the default - see module docstring re: "
                         "point density vs. resonance linewidth")
    p.add_argument("--outdir", required=True)
    p.add_argument("--palace-bin", required=True)
    p.add_argument("--num-cpus", type=int, default=16)
    p.add_argument("--total-length-mm", type=float, default=None,
                    help="Override TOTAL_LENGTH_MM from build_design.py, once "
                         "calibration is finalized this should match whatever "
                         "was used for the eigenmode baseline check.")
    p.add_argument("--bg-min-size-um", type=float, default=15.0,
                    help="Background mesh min_size (um) for TL/resonator1/LP1/LP2, "
                         "away from the grounding strap. Notebook mesh-only testing "
                         "showed 8um->126.6M elements, 15um->19.7M elements (84%% "
                         "reduction) with fillet_resolution making no measurable "
                         "difference either way - this is the dominant lever.")
    p.add_argument("--bg-max-size-um", type=float, default=150.0,
                    help="Background mesh max_size (um), bulk substrate/air-box "
                         "coarsening far from any conductor edge.")
    p.add_argument("--strap-margin-um", type=float, default=3.0,
                    help="Fixed padding (um) around the strap's bounding box for "
                         "the fine-mesh region. Must stay above 0 - a clearance-"
                         "derived version of this previously collapsed to exactly "
                         "0 (since the strap touches its own trace by design) and "
                         "produced a degenerate mesh that crashed Palace instantly.")
    p.add_argument("--solver-order", type=int, default=1,
                    help="FEM element order. order=2 gives better accuracy but "
                         "roughly 4-8x more DOF for the same mesh - an earlier "
                         "test at order=2 with this domain size (feedline + "
                         "both launchpads + resonator) hit an out-of-memory "
                         "kill at ~42M unknowns. order=1 is the safer default "
                         "until node memory limits are better characterized.")
    p.add_argument("--resonator-min-size-um", type=float, default=2.0,
                    help="min mesh element size (um) applied ONLY in a small "
                         "box around the actual grounding strap (extracted "
                         "from the built design's geometry, not the whole "
                         "resonator component). 2um gives ~5 elements across "
                         "the strap's 10um width; the rest of the resonator "
                         "stays at the coarser 8um setting shared with "
                         "TL/LP1/LP2.")
    args = p.parse_args()

    ground_pos_label = "none" if args.ground_pos is None else f"{args.ground_pos}mm"
    total_length_mm = args.total_length_mm if args.total_length_mm is not None else TOTAL_LENGTH_MM
    t0 = time.time()

    freq_step = (FREQ_END_HZ - FREQ_START_HZ) / max(args.n_points - 1, 1)
    print(f"[1/5] Continuous sweep: {FREQ_START_HZ/1e9:.4f} - {FREQ_END_HZ/1e9:.4f} GHz "
          f"({args.n_points} points, step={freq_step/1e6:.2f} MHz)")

    print(f"[2/5] Building design (ground_pos={ground_pos_label}, "
          f"total_length={total_length_mm}mm)...")
    design = build_design(ground_pos_mm=args.ground_pos, total_length_mm=total_length_mm)
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
        "solns_to_save": args.n_points,
        "solver_order": args.solver_order,
        "solver_tol": 1.0e-8,
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

        # Mesh strategy: resonator1 stays at a coarse 8um everywhere (same
        # as TL/LP1/LP2, all well-tested at this scale) EXCEPT for a small
        # box around the actual rendered 10um grounding strap, extracted
        # directly from the built design's qgeometry 'poly' table so this
        # is correct for whatever ground_pos_mm was used - not a guess.
        #
        # This replaces an earlier approach that fine-meshed resonator1's
        # entire ~11.8mm length at 2-4um, which caused two separate OOM
        # kills (55M and 42M unknowns) despite the strap itself only being
        # a ~24x17um feature - the API used (fine_mesh_components) is
        # component-scoped, not feature-scoped, so it was refining the
        # whole trace, not just the strap.
        driven_sim.fine_mesh_components(
            ["TL", "resonator1", "LP1", "LP2"],
            min_size=args.bg_min_size_um * 1e-6, max_size=args.bg_max_size_um * 1e-6,
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

                # Fixed margin (NOT clearance-derived - see below for why).
                # An earlier version measured distance from the strap to the
                # resonator's own trace geometry and used a fraction of that
                # as the margin. This was conceptually broken: the strap is
                # PHYSICALLY ATTACHED to the trace at its position, so that
                # distance is naturally ~0 right there (it's supposed to
                # touch), not a measure of clearance to a genuinely different
                # neighboring meander strand. When that measured distance
                # came out as exactly 0.0, the margin collapsed to 0, GMSH
                # produced a degenerate zero-width fine-mesh region, and
                # Palace crashed almost immediately reading the resulting
                # mesh (signal 9, before any normal mesh-partitioning output
                # even printed - a mesh-validity crash, not an OOM).
                # A fixed, always-positive margin avoids this failure mode
                # entirely. This is also what the one successful run
                # actually used in practice (its clearance-derived value
                # happened to be larger than this fixed cap, so the cap was
                # what got applied anyway).
                safe_margin_mm = args.strap_margin_um / 1000

                x1, y1 = (minx - safe_margin_mm) * 1e-3, (miny - safe_margin_mm) * 1e-3
                x2, y2 = (maxx + safe_margin_mm) * 1e-3, (maxy + safe_margin_mm) * 1e-3
                print(f"      Strap bounding box (mm): x=[{minx:.4f},{maxx:.4f}] "
                      f"y=[{miny:.4f},{maxy:.4f}], margin={safe_margin_mm*1000:.1f}um - "
                      f"applying {args.resonator_min_size_um}um fine mesh there only.")
                driven_sim.fine_mesh_in_rectangle(
                    x1, y1, x2, y2,
                    min_size=args.resonator_min_size_um * 1e-6, max_size=args.bg_min_size_um * 1e-6,
                    taper_dist_min=2e-6, taper_dist_max=6e-6,
                )

        driven_sim.set_freq_values(
            freq_start=FREQ_START_HZ,
            freq_end=FREQ_END_HZ,
            freq_step=freq_step,
        )

        driven_sim.prepare_simulation()
        driven_sim.run()
    except Exception:
        print("\n!!! FAILED during driven simulation !!!")
        traceback.print_exc()
        sys.exit(1)

    print(f"      OK ({time.time()-t1:.1f}s)")

    print("[5/5] Retrieving S-parameter data...")
    s_data = driven_sim.retrieve_data()
    if s_data is None or "freq" not in s_data:
        print("\n!!! FAILED: run() completed without raising, but produced no "
              "usable output (retrieve_data() returned None or incomplete data).")
        print("This usually means the MPI job itself was killed mid-run "
              "(check for OOM/signal 9 above) without Palace/SQDMetal "
              "propagating that as a Python exception.")
        sys.exit(1)

    import numpy as np
    freqs = np.asarray(s_data["freq"])
    s11 = np.asarray(s_data["S11"])
    s21 = np.asarray(s_data["S21"])
    print(f"      OK. Retrieved {len(freqs)} frequency points.")

    out_npz = f"{args.outdir}/s_sweep.npz"
    np.savez(out_npz,
             ground_pos_mm=(-1.0 if args.ground_pos is None else args.ground_pos),
             ground_pos_is_baseline=(args.ground_pos is None),
             freq_hz=freqs,
             S11=s11,
             S21=s21)

    result = {
        "status": "success",
        "ground_pos_mm": args.ground_pos,
        "total_length_mm": total_length_mm,
        "n_points": len(freqs),
        "freq_range_ghz": [FREQ_START_HZ / 1e9, FREQ_END_HZ / 1e9],
        "wall_time_s": time.time() - t0,
    }
    json.dump(result, open(f"{args.outdir}/sweep_summary.json", "w"), indent=2)

    print(f"\n=== STAGE B (pos={ground_pos_label}) COMPLETE in {result['wall_time_s']:.1f}s ===")
    print(f"S-parameter data written to {out_npz}")


if __name__ == "__main__":
    main()
