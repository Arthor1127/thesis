"""Eigenmode check for one (coupling_space_um, coupling_length_um) hanger point,
following this repo's build_and_eigenmode.py convention (see e.g.
scripts/python/transmission_test/build_and_eigenmode.py).

Unlike run_qc_sim.py's driven S21 sim, this has NO ports at all: LP1/LP2 are
left as plain open-ended PEC pads (no create_port_CPW_on_Launcher, no
set_port_excitation), so the launcher-to-launcher feedline itself becomes a
closed structure with its own open-open standing-wave spectrum. That matters
here specifically -- the driven cs14_cl300_ref pilot showed TWO weak resonant
features (7.72GHz and 7.86GHz, see the phase-residual analysis), and one
candidate explanation is that one of them is not the hanger resonator at all
but a bare feedline-cavity mode (LP1<->LP2, ~3.5mm of trace + 2 launch tapers)
that only exists because the driven sim's lumped ports don't behave as a true
open circuit off-resonance. Eigenmode has no such ambiguity: whatever modes
come back in the search window are the geometry's real ones, so this is the
independent cross-check for whether the second driven-sim feature is real
resonator physics or a port artifact.

Cluster job only -- run via a qsub script (this repo's standing rule: no
PALACE/FEM solves locally, mesh-only inspection is fine locally but this is a
real eigensolve). Not wired into a submit_*.sh yet -- run manually:

    python build_and_eigenmode.py --coupling-space-um 14 --coupling-length-um 300 \\
        --outdir out_eig --palace-bin $PALACE_BIN --num-cpus 24 \\
        --amr-max-its 3 --amr-tol 1e-2 --number-of-freqs 4 --starting-freq 6e9
"""
import argparse
import json
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root, for readout_resonators etc.
from readout_resonators import build_design, CPW_WIDTH_UM, CPW_GAP_UM


def fine_mesh_cpw_gap(sim, design, comp_name, trace_name="", dist_resolution_um=20.0,
                      width_um=None, gap_um=None, gap_min_size_um=1.0, gap_max_size_um=80.0,
                      taper_dist_min_um=None, taper_dist_max_um=None):
    """Fine-mesh a thin strip centered on EACH of a CPW path's two lateral
    gaps (offset +-(width/2 + gap/2) from the trace centerline) instead of
    the component's whole footprint. There's no public SQDMetal wrapper for
    this: fine_mesh_along_path only meshes the centerline itself, which for a
    CPW is the metal interior -- the LEAST field-concentrated place, not the
    gap edges where the E-field actually piles up.

    Uses QUtilities.calc_points_on_path (the same helper
    PALACE_Model_Base.fine_mesh_along_path calls internally), which already
    returns the sampled centerline points AND the per-point right-hand normal
    vector plus the trace's own width/gap -- exactly what's needed to build
    the two offset polylines. They're appended directly to sim._fine_meshes
    as 'path' entries, the same list every public fine_mesh_* method
    populates; there's just no method that offsets first.
    """
    from SQDMetal.Utilities.QUtilities import QUtilities
    leUnits = QUtilities.get_units(design)
    pts, normals, path_width, path_gap, _ = QUtilities.calc_points_on_path(
        dist_resolution_um * 1e-6 / leUnits, design, comp_name, trace_name)
    w = width_um * 1e-6 / leUnits if width_um is not None else path_width
    g = gap_um * 1e-6 / leUnits if gap_um is not None else path_gap
    offset = w / 2 + g / 2
    taper_min = (taper_dist_min_um * 1e-6 if taper_dist_min_um is not None
                else gap_min_size_um * 2e-6)
    taper_max = (taper_dist_max_um * 1e-6 if taper_dist_max_um is not None
                else gap_max_size_um * 1e-6)
    for sign in (+1, -1):
        offset_pts = (pts + sign * offset * normals) * leUnits
        sim._fine_meshes.append({
            "type": "path",
            "path": offset_pts,
            "min_size": gap_min_size_um * 1e-6,
            "max_size": gap_max_size_um * 1e-6,
            "taper_dist_min": taper_min,
            "taper_dist_max": taper_max,
        })


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--coupling-space-um", type=float, required=True)
    p.add_argument("--coupling-length-um", type=float, required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--palace-bin", required=True)
    p.add_argument("--num-cpus", type=int, default=16)
    p.add_argument("--starting-freq", type=float, default=6e9)
    p.add_argument("--number-of-freqs", type=int, default=4)
    p.add_argument("--amr-max-its", type=int, default=0)
    p.add_argument("--amr-tol", type=float, default=1e-2)
    p.add_argument("--bg-min-size-um", type=float, default=16.0,
                    help="fine mesh on TL1/TL2/LP1/LP2 (feedline+launchers). Doubled from "
                        "8um: the feedline is now 12mm long (vs 2mm originally), so this "
                        "resolution applies along 6x the length -- the local mesh-only "
                        "check came back a big 59MB even before any cluster solve.")
    p.add_argument("--bg-max-size-um", type=float, default=80.0,
                    help="air/dielectric background max element size. Lowered from an "
                        "earlier 150um default: with launchpads at the default 4mm-chip's "
                        "+-1mm and only ~100um margin to the chip edge, and the fine "
                        "component mesh going down to ~1-7um, 150um left too coarse a "
                        "transition zone between component mesh and background -- see "
                        "--chip-size-x-mm/--chip-size-y-mm below for the matching geometry fix.")
    p.add_argument("--resonator-min-size-um", type=float, default=None,
                    help="fine mesh on readout_1/CLT/res_short. "
                        "Defaults to 0.5*min(cpw_width,cpw_gap,coupling_space) "
                        "same rule run_qc_sim.py uses, floor 1um.")
    p.add_argument("--chip-half-span-um", type=float, default=None,
                    help="LP1/LP2 distance from centre = half the launcher-to-launcher "
                        "feedline length. Defaults to readout_resonators.CHIP_HALF_SPAN_UM. "
                        "This is the swept variable in submit_feedline_sweep.sh: bare "
                        "feedline harmonics scale as v/(2L) so they move with it, while the "
                        "resonator's own mode does not -- which is what tells the two apart.")
    p.add_argument("--chip-size-x-mm", type=float, default=None,
                    help="Defaults to None -> build_design auto-sizes it to the launchpads "
                        "plus margin, which is what a feedline sweep wants (the chip has to "
                        "grow with the feedline).")
    p.add_argument("--chip-size-y-mm", type=float, default=4.0,
                    help="Matches build_design's own current default. The resonator is now "
                        "a single east-jog-west fold (SNAKE_X_MAX_UM raised way up since the "
                        "launchpads are much farther away) reaching down to at most "
                        "~-1.2mm (worst case: coupling_length_um=100, down_length_um=600); "
                        "4mm (+-2mm) keeps ~800um of margin past that.")
    args = p.parse_args()

    from readout_resonators import CHIP_HALF_SPAN_UM
    half_span_um = (args.chip_half_span_um if args.chip_half_span_um is not None
                    else CHIP_HALF_SPAN_UM)

    os.makedirs(args.outdir, exist_ok=True)
    t0 = time.time()
    print(f"[1/4] Building design (coupling_space={args.coupling_space_um}um, "
          f"coupling_length={args.coupling_length_um}um, "
          f"half_span={half_span_um}um -> feedline={2*half_span_um/1000:.1f}mm)...", flush=True)
    design = build_design(args.coupling_space_um, args.coupling_length_um,
                          chip_size_x_mm=args.chip_size_x_mm, chip_size_y_mm=args.chip_size_y_mm,
                          chip_half_span_um=half_span_um)
    design.rebuild()
    print(f"      chip = {design.chips.main.size.size_x} x {design.chips.main.size.size_y}",
          flush=True)
    print(f"      OK ({time.time()-t0:.1f}s). Components: {list(design.components.keys())}",
          flush=True)

    print("[2/4] Importing SQDMetal PALACE eigenmode module...", flush=True)
    from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
    import matplotlib
    matplotlib.use("Agg", force=True)
    print("      OK", flush=True)

    resonator_min_um = (args.resonator_min_size_um if args.resonator_min_size_um is not None
                        else max(0.5 * min(CPW_WIDTH_UM, CPW_GAP_UM, args.coupling_space_um), 1.0))

    print("[3/4] Setting up + running eigenmode simulation...", flush=True)
    t1 = time.time()
    user_defined_options = {
        "mesh_refinement": 0,  # UniformLevels - NOT AMR; see --amr-max-its
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
            name=f"eig_cs{args.coupling_space_um:g}_cl{args.coupling_length_um:g}"
                 f"_hs{half_span_um:g}",
            metal_design=design,
            sim_parent_directory=args.outdir,
            mode="simPC",
            meshing="GMSH",
            user_options=user_defined_options,
            create_files=True,
        )
        eigen_sim.add_metallic(1)
        eigen_sim.add_ground_plane()
        # NO ports -- LP1/LP2 stay plain open-ended PEC pads on purpose (see
        # module docstring): this lets the bare feedline's own open-open
        # cavity modes show up alongside the hanger resonator's mode, which
        # is exactly the ambiguity we're trying to resolve.
        eigen_sim.fine_mesh_components(
            ["TL1", "TL2", "LP1", "LP2"], min_size=args.bg_min_size_um * 1e-6,
            max_size=args.bg_max_size_um * 1e-6, taper_dist_min=16e-6, taper_dist_max=40e-6,
            metals_only=False,
        )
        # readout_1 is ~3.7mm of routed trace -- fine-meshing its WHOLE footprint
        # (metals_only=False, as the old fine_mesh_components call below still
        # does for CLT/res_short) burns a huge number of elements on plain
        # meander interior that has no field concentration at all. Only the two
        # thin strips running along its own CPW gap actually need resonator_min_um
        # resolution; fine_mesh_cpw_gap targets exactly those (see its docstring).
        fine_mesh_cpw_gap(
            eigen_sim, design, "readout_1", dist_resolution_um=20.0,
            gap_min_size_um=resonator_min_um, gap_max_size_um=args.bg_max_size_um,
            taper_dist_min_um=resonator_min_um * 2, taper_dist_max_um=resonator_min_um * 4,
        )
        # CLT and res_short are tiny (order coupling_length_um + termination stub) --
        # full-footprint fine meshing there is cheap, so left as-is rather than
        # guessing at CoupledLineTee's internal coupled-line gap geometry.
        eigen_sim.fine_mesh_components(
            ["CLT", "res_short"], min_size=resonator_min_um * 1e-6,
            max_size=args.bg_max_size_um * 1e-6, taper_dist_min=30e-6, taper_dist_max=100e-6,
            metals_only=False,
        )
        print(f"      mesh: bg_min={args.bg_min_size_um}um resonator_min={resonator_min_um:.3f}um",
              flush=True)
        # setup_EPR_interfaces() deliberately omitted -- it exists to support energy
        # participation ratio / interface-loss postprocessing, which we don't need
        # here (we only want mode frequencies + Qs). It also measurably raises mesh
        # complexity: job 1938821's log showed it adding 19590 duplicate vertices,
        # 44223 duplicate boundary elements, and 902 material-interface boundary
        # elements from splitting the mesh at every metal-air/substrate-air/
        # substrate-metal interface -- extra DOFs and extra degenerate-element risk
        # for physics this script doesn't use.
        if args.amr_max_its > 0:
            eigen_sim._mesh_refinement = {
                "UniformLevels": 0,
                "MaxIts": args.amr_max_its,
                "Tol": args.amr_tol,
            }
            print(f"      Solution-based AMR ENABLED: MaxIts={args.amr_max_its}, "
                  f"Tol={args.amr_tol:.1e}", flush=True)

        eigen_sim.prepare_simulation()
        eigen_sim.run()
    except Exception:
        print("\n!!! FAILED during eigenmode simulation !!!", flush=True)
        traceback.print_exc()
        json.dump({"status": "failed", "coupling_space_um": args.coupling_space_um,
                   "coupling_length_um": args.coupling_length_um,
                   "chip_half_span_um": half_span_um},
                   open(f"{args.outdir}/resonance.json", "w"), indent=2)
        sys.exit(1)

    print(f"      OK ({time.time()-t1:.1f}s)", flush=True)

    print("[4/4] Parsing eig.csv output...", flush=True)
    eig_csv_path = eigen_sim._output_data_dir + "/eig.csv"
    if not os.path.exists(eig_csv_path):
        print(f"\n!!! FAILED: run() completed without raising, but {eig_csv_path} "
              "was never written (likely an MPI OOM-kill mid-run).", flush=True)
        json.dump({"status": "failed", "coupling_space_um": args.coupling_space_um,
                   "coupling_length_um": args.coupling_length_um,
                   "reason": "missing eig.csv, likely killed MPI job"},
                   open(f"{args.outdir}/resonance.json", "w"), indent=2)
        sys.exit(1)

    import pandas as pd
    eig_df = pd.read_csv(eig_csv_path)
    eig_df.columns = eig_df.columns.str.strip()
    f_re_col = [c for c in eig_df.columns if c.startswith("Re{f}")][0]
    q_col = [c for c in eig_df.columns if c.startswith("Q")][0]

    if len(eig_df) == 0:
        print(f"\n!!! FAILED: eig.csv has a header but ZERO converged modes in the "
              f"search window (starting_freq={args.starting_freq/1e9:.3f}GHz, "
              f"number_of_freqs={args.number_of_freqs}).", flush=True)
        json.dump({"status": "no_modes_converged", "coupling_space_um": args.coupling_space_um,
                   "coupling_length_um": args.coupling_length_um,
                   "starting_freq_hz": args.starting_freq,
                   "number_of_freqs": args.number_of_freqs},
                   open(f"{args.outdir}/resonance.json", "w"), indent=2)
        sys.exit(1)

    all_modes_ghz = eig_df[f_re_col].astype(float).tolist()
    all_qs = eig_df[q_col].astype(float).tolist()
    print(f"      All modes found (GHz): {all_modes_ghz}", flush=True)
    print(f"      Q factors:            {all_qs}", flush=True)

    result = {
        "status": "success",
        "coupling_space_um": args.coupling_space_um,
        "coupling_length_um": args.coupling_length_um,
        "chip_half_span_um": half_span_um,
        "feedline_length_mm": 2 * half_span_um / 1000.0,
        "all_modes_ghz": all_modes_ghz,
        "all_q_factors": all_qs,
        "starting_freq_hz": args.starting_freq,
        "number_of_freqs": args.number_of_freqs,
        "amr_max_its": args.amr_max_its,
        "amr_tol": args.amr_tol,
        "num_cpus": args.num_cpus,
        "wall_time_s": time.time() - t0,
    }
    out_path = f"{args.outdir}/resonance.json"
    json.dump(result, open(out_path, "w"), indent=2)
    print(f"\n=== DONE in {result['wall_time_s']:.1f}s === Result written to {out_path}",
          flush=True)


if __name__ == "__main__":
    main()
