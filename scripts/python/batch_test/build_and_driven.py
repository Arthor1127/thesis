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

--- On mesh background sizing ---
Ports (TL/LP1/LP2) and resonator1 have INDEPENDENT background mesh sizes
(--bg-min/max-size-um for ports, --resonator-bg-min/max-size-um for the
resonator, which defaults to the port size if omitted). resonator1 is by
far the largest component by trace length, so this split is the main
lever for cutting element/DOF count without touching port resolution or
--solver-order. Example, coarsening just the resonator:
    python build_and_driven.py --ground-pos 5.906 \
        --outdir ~/sweep_v2/driven_5.906 --palace-bin $PALACE_BIN --num-cpus 16 \
        --bg-min-size-um 15 --bg-max-size-um 150 \
        --resonator-bg-min-size-um 30 --resonator-bg-max-size-um 300

Usage - validated production settings (order=2, split mesh, sane
solns/tol) - this is what stageB_driven.sh actually runs:
    python build_and_driven.py --ground-pos 5.906 \
        --outdir ~/sweep_v2/driven_5.906 --palace-bin $PALACE_BIN --num-cpus 16 \
        --bg-min-size-um 15 --bg-max-size-um 150 \
        --resonator-bg-min-size-um 30 --resonator-bg-max-size-um 300 \
        --solver-order 2 --solns-to-save 10 --solver-tol 1e-8
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
                    help="Background mesh min_size (um) for TL/LP1/LP2 (the ports), "
                         "which carry the actual 50 ohm lumped-port excitation and "
                         "matching, so they keep the finer/validated resolution. "
                         "Notebook mesh-only testing on the original (shorter, "
                         "port+resonator combined) mesh showed 8um->126.6M elements, "
                         "15um->19.7M elements (84%% reduction) - background size is "
                         "the dominant memory lever.")
    p.add_argument("--bg-max-size-um", type=float, default=150.0,
                    help="Background mesh max_size (um) for TL/LP1/LP2, bulk "
                         "substrate/air-box coarsening far from any conductor edge.")
    p.add_argument("--resonator-bg-min-size-um", type=float, default=None,
                    help="Background mesh min_size (um) for resonator1 specifically, "
                         "independent of --bg-min-size-um. resonator1 is a long "
                         "(mm-scale) CPW trace propagating a mode already understood "
                         "analytically away from the strap; unlike TL/LP1/LP2 it does "
                         "not need port-grade resolution along its whole length, only "
                         "at the strap (handled separately via fine_mesh_in_rectangle "
                         "below). Splitting this out from --bg-min-size-um lets the "
                         "resonator go much coarser without touching validated port "
                         "behavior. Defaults to --bg-min-size-um if not given (i.e. "
                         "old single-value behavior).")
    p.add_argument("--resonator-bg-max-size-um", type=float, default=None,
                    help="Background mesh max_size (um) for resonator1 specifically. "
                         "Defaults to --bg-max-size-um if not given.")
    p.add_argument("--chip-size-x-mm", type=float, default=None,
                    help="Chip X footprint (mm). None uses build_design.py's "
                         "default. Little slack here: launchpads reach ~+-2.1mm "
                         "against a 4.8mm chip's +-2.4mm edge.")
    p.add_argument("--chip-size-y-mm", type=float, default=None,
                    help="Chip Y footprint (mm). None uses build_design.py's "
                         "default (2.4). This is where the real slack is - "
                         "~0.8mm below the lowest feature is dead substrate. "
                         "1.8 (with --chip-center-y-mm -0.7) keeps 250um of "
                         "ground below otg2 for ~75%% of the volume. Element "
                         "count, and so AMG setup time, scales with volume. "
                         "CHANGING THIS SHIFTS f1 - re-run the ungrounded "
                         "baseline before trusting a sweep.")
    p.add_argument("--chip-center-x-mm", type=float, default=None,
                    help="Chip X center (mm). None uses the module default.")
    p.add_argument("--chip-center-y-mm", type=float, default=None,
                    help="Chip Y center (mm). None uses the module default. "
                         "Pair with --chip-size-y-mm to keep the content "
                         "centred while shrinking (e.g. 1.8 / -0.7).")
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
                         "until node memory limits are better characterized. "
                         "NOTE: production runs so far have used order=2 "
                         "explicitly (--solver-order 2) on large-memory nodes - "
                         "stageB_driven.sh passes this explicitly, do not rely "
                         "on this default silently changing that.")
    p.add_argument("--adaptive-tol", type=float, default=1e-3,
                    help="Relative error tolerance for Palace's adaptive fast "
                         "frequency sweep (PROM). >0 enables it: Palace solves "
                         "the full-order model at only a few auto-selected "
                         "sample frequencies, builds a reduced-order model, and "
                         "evaluates all remaining points from that ROM cheaply. "
                         "Set to 0 to disable and force a full-order solve at "
                         "every frequency (the old, slow behavior). NOTE: this "
                         "is an L2 coefficient-space indicator, not a strict "
                         "error bound on S-parameters - validate against a few "
                         "uniform points on a new design before trusting it.")
    p.add_argument("--adaptive-max-samples", type=int, default=30,
                    help="Cap on how many full-order solves the adaptive sweep "
                         "may use if --adaptive-tol is not reached first. "
                         "Palace's own default is 20; 30 gives headroom for a "
                         "wide 4.5-20.5GHz band spanning four harmonics, since "
                         "each resonance needs samples near it to be captured "
                         "by the ROM.")
    p.add_argument("--solns-to-save", type=int, default=10,
                    help="Number of frequency points to write full 3D field "
                         "snapshots for (disk I/O only, does not affect solve "
                         "memory/time). Defaults to a small number, NOT "
                         "--n-points - writing a full field snapshot at all "
                         "321 points x 10 sweep positions is unnecessary disk "
                         "churn when you only need a handful for visualization.")
    p.add_argument("--solver-tol", type=float, default=1.0e-8,
                    help="FGMRES convergence tolerance. 1e-8 is tight - a real "
                         "production point took 106 iterations to converge at "
                         "this tolerance. 1e-6 is very likely sufficient for "
                         "S-parameter engineering accuracy and should reduce "
                         "iterations-per-point meaningfully; not yet validated "
                         "against 1e-8 results on this design, so treat this as "
                         "an experiment to check, not an assumed-safe default.")
    p.add_argument("--resonator-min-size-um", type=float, default=2.0,
                    help="min mesh element size (um) applied ONLY in a small "
                         "box around the actual grounding strap (extracted "
                         "from the built design's geometry, not the whole "
                         "resonator component). 2um gives ~5 elements across "
                         "the strap's 10um width; the rest of the resonator "
                         "stays at the (now independently-set) resonator "
                         "background size - see --resonator-bg-min-size-um.")
    args = p.parse_args()

    # resonator1's background defaults to the port background if not given
    # separately, preserving old single-value behavior when these flags
    # are omitted.
    resonator_bg_min_um = (args.resonator_bg_min_size_um
                            if args.resonator_bg_min_size_um is not None
                            else args.bg_min_size_um)
    resonator_bg_max_um = (args.resonator_bg_max_size_um
                            if args.resonator_bg_max_size_um is not None
                            else args.bg_max_size_um)

    ground_pos_label = "none" if args.ground_pos is None else f"{args.ground_pos}mm"
    total_length_mm = args.total_length_mm if args.total_length_mm is not None else TOTAL_LENGTH_MM
    t0 = time.time()

    freq_step = (FREQ_END_HZ - FREQ_START_HZ) / max(args.n_points - 1, 1)
    print(f"[1/5] Continuous sweep: {FREQ_START_HZ/1e9:.4f} - {FREQ_END_HZ/1e9:.4f} GHz "
          f"({args.n_points} points, step={freq_step/1e6:.2f} MHz)")
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

        # Mesh strategy, split by role:
        #
        # - TL/LP1/LP2 carry the actual 50 ohm lumped-port excitation and
        #   impedance matching, so they keep the finer/validated background
        #   size (--bg-min/max-size-um).
        # - resonator1 is a long (mm-scale) CPW trace whose physics away
        #   from the strap is already understood analytically - it does
        #   NOT need port-grade resolution along its whole length, only at
        #   the strap itself (handled separately below). It gets its own,
        #   independently coarser background size (--resonator-bg-min/max
        #   -size-um, defaults to the port size if not overridden).
        # - A small box around the actual rendered 10um grounding strap,
        #   extracted directly from the built design's qgeometry 'poly'
        #   table (correct for whatever ground_pos_mm was used, not a
        #   guess), gets the finest resolution of all.
        #
        # This two-call split replaces the original single fine_mesh_
        # components() call over ["TL","resonator1","LP1","LP2"] together,
        # which forced resonator1 (by far the largest component by trace
        # length) to share the ports' resolution across its entire length
        # - the dominant driver of element/DOF count on this design. An
        # earlier, cruder attempt to avoid this fine-meshed resonator1's
        # entire ~11.8mm length at 2-4um directly (instead of splitting
        # backgrounds), which caused two separate OOM kills (55M and 42M
        # unknowns) - fine_mesh_components() is component-scoped, not
        # feature-scoped, so it refines the whole trace it's given, not
        # just a feature on it.
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
                # max_size here is the taper TARGET - what the strap's fine
                # mesh transitions back out to. Now that resonator1 has its
                # own (independently coarser) background, that target must
                # be resonator_bg_min_um, not the ports' bg_min_size_um -
                # otherwise the strap box would taper toward a resolution
                # finer than resonator1's actual surrounding mesh, forcing
                # an unnecessary extra refinement step right at the box
                # boundary.
                driven_sim.fine_mesh_in_rectangle(
                    x1, y1, x2, y2,
                                    # max_size here is the taper TARGET - the size this
                # field relaxes to far from the box. It MUST be the global
                # background MAX, not the min. GMSH's threshold field
                # returns SizeMax as a VALUE at every point beyond
                # taper_dist_max (it does not mean "unconstrained there"),
                # and GMSH takes the pointwise MINIMUM over all size
                # fields - so setting this to the background MIN silently
                # caps the ENTIRE domain at that size and destroys the
                # far-field coarsening. Measured on the quarter-wave
                # design: baseline 60,110 elements vs 7,810,000 with a
                # strap at ANY position (+-0.14% across 6 positions, i.e.
                # completely position-independent), matching a uniform
                # background-min mesh over the whole domain.
                    min_size=args.resonator_min_size_um * 1e-6, max_size=resonator_bg_max_um * 1e-6,
                    taper_dist_min=2e-6, taper_dist_max=6e-6,
                )

        driven_sim.set_freq_values(
            freq_start=FREQ_START_HZ,
            freq_end=FREQ_END_HZ,
            freq_step=freq_step,
        )

        driven_sim.prepare_simulation()

        # --- Adaptive fast frequency sweep (PROM) ---
        # By default Palace solves the FULL-ORDER model at every one of the
        # n_points frequencies - 321 complete solves per position, each
        # rebuilding the AMG preconditioner, which is ~75% of runtime. With
        # AdaptiveTol > 0 Palace instead solves the full-order model at a
        # few automatically-chosen sample frequencies (capped by
        # AdaptiveMaxSamples, default 20), builds a projection-based
        # reduced-order model from them, and evaluates the remaining points
        # from that ROM at negligible cost. For a 321-point sweep that is
        # potentially a ~16x reduction in expensive solves.
        #
        # SQDMetal does not expose these keys (its set_freq_values only
        # writes MinFreq/MaxFreq/FreqStep), so patch the config JSON
        # directly here. Same technique SQDMetal itself uses in
        # Frequency_Driven_Simulation.set_freq_values(), and it must happen
        # AFTER prepare_simulation() (which writes the config) but BEFORE
        # run() (which consumes it).
        #
        # CAVEAT worth knowing: per Palace's own reference docs, AdaptiveTol
        # is an L2 coefficient-space indicator, NOT a strict relative-error
        # bound on derived quantities like S-parameters. Validate a new
        # design against a handful of uniformly-sampled points before
        # trusting an adaptive sweep wholesale.
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
    # SQDMetal's PALACE_Driven_Simulation.retrieve_data() returns a dict
    # keyed 'freqs' (plural) - not 'freq'. An earlier version of this
    # script checked/accessed "freq" (singular), which meant a genuinely
    # successful Palace run (real port-S.csv written, retrieve_data()
    # returning a valid dict) could still get reported as FAILED here
    # purely from this key-name mismatch, with no actual data problem.
    # Confirmed by inspecting Frequency_Driven_Simulation.py directly:
    # `ret_data = {'freqs': freq_vals}`.
    if s_data is None or "freqs" not in s_data:
        print("\n!!! FAILED: run() completed without raising, but produced no "
              "usable output (retrieve_data() returned None or incomplete data).")
        print("This usually means the MPI job itself was killed mid-run "
              "(check for OOM/signal 9 above) without Palace/SQDMetal "
              "propagating that as a Python exception.")
        sys.exit(1)

    import numpy as np
    freqs = np.asarray(s_data["freqs"])
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
