import argparse
import json
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import (build_design, N_SQUIDS, L_J_NH, CJ_J_FF, RES_LENGTH_MM,
                          JJ_PER_LOOP, ARRAY_COMPONENT_NAME,
                          SQUID_GAP_UM, ARM_WIDTH_UM, JOSEPHSON_GAP_UM,
                          SEGMENT_LENGTH_UM, ARM_SEMI_LENGTH_UM, LINE_WIDTH_UM,
                          GROUND_GAP_UM, ARRAY_POS_X_UM, ARRAY_POS_Y_UM)


def _json_safe(obj):
    """Recursively convert numpy arrays/scalars into plain Python types so
    json.dump doesn't choke on them.

    numpy complex scalars have .tolist(), but it converts them to a plain
    Python `complex` - which is STILL not JSON serializable, so that result
    must be re-checked (this once crashed mid-write on interface_epr data,
    corrupting the output file for a run whose solve had already succeeded)."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, complex):
        return {"real": obj.real, "imag": obj.imag}
    if hasattr(obj, "tolist"):
        return _json_safe(obj.tolist())
    return obj


def _build_parser():
    p = argparse.ArgumentParser(
        description="Build the SQUID-array + readout-resonator design and run a "
                    "Palace eigenmode EPR simulation. Every CLI knob below that "
                    "defaults to None falls back to the constant of the same "
                    "name in build_design.py.")

    # --- what to launch: the physics sweep knobs ---
    sweep = p.add_argument_group("sweep knobs (geometry / inductance)")
    sweep.add_argument("--n-squids", type=int, default=None,
                       help=f"Number of SQUID loops in the series array. "
                            f"Default build_design.N_SQUIDS ({N_SQUIDS}).")
    sweep.add_argument("--L-j-nH", type=float, default=None,
                       help=f"Inductance of EVERY junction in the array (nH). "
                            f"Default build_design.L_J_NH ({L_J_NH}). Total array "
                            f"superinductance ~ n_squids * L_j / 2 (two junctions "
                            f"per loop in parallel).")
    sweep.add_argument("--Cj-j-fF", type=float, default=None,
                       help=f"Parallel capacitance of every junction (fF). "
                            f"Default build_design.CJ_J_FF ({CJ_J_FF}).")
    sweep.add_argument("--res-length-mm", type=float, default=None,
                       help=f"Readout resonator total length (mm). "
                            f"Default build_design.RES_LENGTH_MM ({RES_LENGTH_MM}).")
    sweep.add_argument("--squid-gap-um", type=float, default=None,
                       help=f"Loop height (um). Default {SQUID_GAP_UM}.")
    sweep.add_argument("--arm-width-um", type=float, default=None,
                       help=f"Width of plate1 / seg a / seg b / seg c (um). "
                            f"Default {ARM_WIDTH_UM}.")
    sweep.add_argument("--josephson-gap-um", type=float, default=None,
                       help=f"Width of the gap each junction bridges (um). "
                            f"Default {JOSEPHSON_GAP_UM}.")
    sweep.add_argument("--segment-length-um", type=float, default=None,
                       help=f"seg d (right contact) length (um). Default {SEGMENT_LENGTH_UM}.")
    sweep.add_argument("--arm-semi-length-um", type=float, default=None,
                       help=f"seg a / seg b length (um). Default {ARM_SEMI_LENGTH_UM}.")
    sweep.add_argument("--line-width-um", type=float, default=None,
                       help=f"seg d (right contact) height (um). Default {LINE_WIDTH_UM}.")
    sweep.add_argument("--ground-gap-um", type=float, default=None,
                       help=f"Ground-plane pocket width around the whole chain (um). "
                            f"Default {GROUND_GAP_UM}.")
    sweep.add_argument("--array-pos-x-um", type=float, default=None,
                       help=f"Array anchor x (um). Default {ARRAY_POS_X_UM}.")
    sweep.add_argument("--array-pos-y-um", type=float, default=None,
                       help=f"Array anchor y (um). Default {ARRAY_POS_Y_UM}.")

    # --- solver ---
    solver = p.add_argument_group("eigenmode solver")
    solver.add_argument("--number-of-freqs", type=int, default=6,
                        help="Number of eigenmodes to solve (resonator fundamental "
                             "+ harmonics). >=6 recommended.")
    solver.add_argument("--starting-freq", type=float, default=4e9,
                        help="Eigen-search start frequency (Hz). The loaded array "
                             "pulls the resonator DOWN, so this is below the "
                             "unloaded lambda/4 frequency.")
    solver.add_argument("--solver-order", type=int, default=1)
    solver.add_argument("--solver-tol", type=float, default=1.0e-5)
    solver.add_argument("--solver-maxits", type=int, default=300)
    solver.add_argument("--fillet-resolution", type=int, default=12)
    solver.add_argument("--amr-max-its", type=int, default=0,
                        help="Solution-based AMR iterations. 0 = off (production).")
    solver.add_argument("--amr-tol", type=float, default=1e-2)

    # --- mesh ---
    mesh = p.add_argument_group("mesh")
    mesh.add_argument("--bg-min-size-um", type=float, default=6.0,
                      help="General background min element size (um); also the "
                           "fallback for the resonator/array min sizes.")
    mesh.add_argument("--bg-max-size-um", type=float, default=150.0)
    mesh.add_argument("--resonator-bg-min-size-um", type=float, default=None,
                      help="readout_1 fine-mesh min size (um). Should resolve "
                           "cpw_gap=12.25um - ~6um is a reasonable floor. "
                           "Defaults to --bg-min-size-um.")
    mesh.add_argument("--resonator-bg-max-size-um", type=float, default=None,
                      help="readout_1 fine-mesh max size (um). Defaults to --bg-max-size-um.")
    mesh.add_argument("--resonator-taper-dist-min-um", type=float, default=40.0,
                      help="readout_1: distance from the metal held at min size (um).")
    mesh.add_argument("--resonator-taper-dist-max-um", type=float, default=200.0,
                      help="readout_1: distance from the metal at which mesh reaches max size (um).")
    mesh.add_argument("--array-bg-min-size-um", type=float, default=2.0,
                      help="Array fine-mesh min size (um). The finest features are "
                           "the josephson gap and ground_gap (~3um each); ~half of "
                           "that resolves them.")
    mesh.add_argument("--array-bg-max-size-um", type=float, default=150.0,
                      help="Array fine-mesh max size (um).")
    mesh.add_argument("--array-taper-dist-min-um", type=float, default=5.0,
                      help="array: distance from the metal held at min size (um). The "
                           "junction-gap current (hence p, hence c ~ sqrt(p)) is set by "
                           "the field in the VACUUM around the 3um gap, so pushing this "
                           "out matters as much as a smaller min size.")
    mesh.add_argument("--array-taper-dist-max-um", type=float, default=200.0,
                      help="array: distance from the metal at which mesh reaches max size (um).")
    mesh.add_argument("--coupling-rect-min-um", type=float, default=None,
                      help="If set, add fine_mesh_in_rectangle boxes around the two "
                           "signal loops (array near end + far end, each with a margin "
                           "of vacuum) at this min size. This is where p - hence c - "
                           "lives; a global fine taper over the whole 5mm of metal "
                           "perimeter just explodes the mesh.")
    mesh.add_argument("--coupling-rect-margin-um", type=float, default=180.0,
                      help="half-size of the coupling-zone boxes (um).")

    # --- chip / infra ---
    chip = p.add_argument_group("chip / infra")
    chip.add_argument("--chip-size-x-mm", type=float, default=None)
    chip.add_argument("--chip-size-y-mm", type=float, default=None)
    chip.add_argument("--chip-center-x-mm", type=float, default=None)
    chip.add_argument("--chip-center-y-mm", type=float, default=None)
    chip.add_argument("--outdir", required=True)
    chip.add_argument("--palace-bin", required=True)
    chip.add_argument("--num-cpus", type=int, default=16)

    return p


def main():
    args = _build_parser().parse_args()

    def _or_default(v, default):
        return default if v is None else v

    n_squids = _or_default(args.n_squids, N_SQUIDS)
    L_j_nH = _or_default(args.L_j_nH, L_J_NH)
    cj_j_fF = _or_default(args.Cj_j_fF, CJ_J_FF)
    res_length_mm = _or_default(args.res_length_mm, RES_LENGTH_MM)

    geom = dict(
        squid_gap_um=_or_default(args.squid_gap_um, SQUID_GAP_UM),
        arm_width_um=_or_default(args.arm_width_um, ARM_WIDTH_UM),
        josephson_gap_um=_or_default(args.josephson_gap_um, JOSEPHSON_GAP_UM),
        segment_length_um=_or_default(args.segment_length_um, SEGMENT_LENGTH_UM),
        arm_semi_length_um=_or_default(args.arm_semi_length_um, ARM_SEMI_LENGTH_UM),
        line_width_um=_or_default(args.line_width_um, LINE_WIDTH_UM),
        ground_gap_um=_or_default(args.ground_gap_um, GROUND_GAP_UM),
        array_pos_x_um=_or_default(args.array_pos_x_um, ARRAY_POS_X_UM),
        array_pos_y_um=_or_default(args.array_pos_y_um, ARRAY_POS_Y_UM),
    )

    os.makedirs(args.outdir, exist_ok=True)
    n_ports = JJ_PER_LOOP * n_squids

    t0 = time.time()
    print(f"[1/4] Building design (n_squids={n_squids}, L_j={L_j_nH}nH, "
          f"res_length={res_length_mm}mm)...")
    design = build_design(n_squids=n_squids, L_j_nH=L_j_nH, cj_j_fF=cj_j_fF,
                          res_length_mm=res_length_mm,
                          chip_size_x_mm=args.chip_size_x_mm,
                          chip_size_y_mm=args.chip_size_y_mm,
                          chip_center_x_mm=args.chip_center_x_mm,
                          chip_center_y_mm=args.chip_center_y_mm,
                          **geom)
    print(f"      OK ({time.time()-t0:.1f}s). Components: {list(design.components.keys())}")

    design.rebuild()

    n_jj = len(design.qgeometry.tables['junction'].loc[
        design.qgeometry.tables['junction']['component']
        == design.components[ARRAY_COMPONENT_NAME].id])
    if n_jj != n_ports:
        print(f"\n!!! FAILED: array has {n_jj} junction rows, expected {n_ports} "
              f"(= {JJ_PER_LOOP} x {n_squids}).")
        json.dump({"status": "failed", "reason": "junction count mismatch",
                   "n_junction_rows": n_jj, "expected": n_ports},
                  open(f"{args.outdir}/epr_result.json", "w"), indent=2)
        sys.exit(1)

    print("[2/4] Importing SQDMetal PALACE eigenmode module...")
    from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
    from SQDMetal.Utilities.Materials import MaterialInterface
    print("      OK")

    import matplotlib
    matplotlib.use("Agg", force=True)

    print(f"[3/4] Setting up + running eigenmode simulation "
          f"({n_ports} lumped inductive ports on the array)...")
    t1 = time.time()
    user_defined_options = {
        "mesh_refinement": 0,   # UniformLevels - NOT AMR; see --amr-max-its
        "dielectric_material": "silicon",
        "starting_freq": args.starting_freq,
        "number_of_freqs": args.number_of_freqs,
        "solns_to_save": args.number_of_freqs,
        "solver_order": args.solver_order,
        "solver_tol": args.solver_tol,
        "solver_maxits": args.solver_maxits,
        "fillet_resolution": args.fillet_resolution,
        "palace_dir": args.palace_bin,
        "num_cpus": args.num_cpus,
    }

    run_label = f"nsq_{n_squids}_Lj_{L_j_nH}nH_reslen_{res_length_mm}mm"

    def _fail(status, extra=None):
        payload = {"status": status, "n_squids": n_squids, "l_j_nh": L_j_nH,
                   "res_length_mm": res_length_mm}
        if extra:
            payload.update(extra)
        json.dump(payload, open(f"{args.outdir}/epr_result.json", "w"), indent=2)

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

        # One lumped inductive port per junction gap. SquidArray emits, per
        # loop i: jj_i_upper (index 2i) then jj_i_lower (index 2i+1). Every gap
        # is a real vacuum break and MUST be ported - an open gap breaks the
        # loop and its participation is identically zero
        # (../squid_resonator_epr/EPR_C_EXTRACTION.md Sec 2). port-EPR.csv will
        # then have n_ports participation columns, in this same index order.
        for k in range(n_ports):
            eigen_sim.create_port_JosephsonJunction(
                ARRAY_COMPONENT_NAME, junction_index=k,
                L_J=L_j_nH * 1e-9, C_J=cj_j_fF * 1e-15)

        resonator_bg_min_um = _or_default(args.resonator_bg_min_size_um, args.bg_min_size_um)
        resonator_bg_max_um = _or_default(args.resonator_bg_max_size_um, args.bg_max_size_um)

        eigen_sim.fine_mesh_components(
            ["readout_1"],
            min_size=resonator_bg_min_um * 1e-6, max_size=resonator_bg_max_um * 1e-6,
            taper_dist_min=args.resonator_taper_dist_min_um * 1e-6,
            taper_dist_max=args.resonator_taper_dist_max_um * 1e-6, metals_only=False,
        )
        eigen_sim.fine_mesh_components(
            [ARRAY_COMPONENT_NAME],
            min_size=args.array_bg_min_size_um * 1e-6, max_size=args.array_bg_max_size_um * 1e-6,
            taper_dist_min=args.array_taper_dist_min_um * 1e-6,
            taper_dist_max=args.array_taper_dist_max_um * 1e-6, metals_only=False,
        )

        if args.coupling_rect_min_um is not None:
            # the array (SquidArray) is anchored at pos_x = ARRAY_POS_X_UM and
            # runs +x; its loop pitch is ~ full_w. The two loops that carry the
            # coupling are the near end (last loop, ~ pin b) and the far end
            # (first loop). Box each with a vacuum margin. (ARRAY_POS_*_UM are
            # imported at module top.)
            pitch_um = 47.0   # SquidArray full_w for the default single-loop geometry
            x_far = ARRAY_POS_X_UM
            x_near = ARRAY_POS_X_UM + n_squids * pitch_um
            m = args.coupling_rect_margin_um
            for tag, xc in (("far", x_far), ("near", x_near)):
                eigen_sim.fine_mesh_in_rectangle(
                    (xc - m) * 1e-6, (ARRAY_POS_Y_UM - m) * 1e-6,
                    (xc + m) * 1e-6, (ARRAY_POS_Y_UM + m) * 1e-6,
                    mesh_min=args.coupling_rect_min_um * 1e-6, mesh_max=60e-6,
                    taper_dist_min=40e-6, taper_dist_max=150e-6)
            print(f"      coupling-zone fine mesh: {args.coupling_rect_min_um}um boxes "
                  f"(+/-{m}um) at x={x_far:.0f} and x={x_near:.0f} um")

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

    print("[4/4] Retrieving EPR data and writing epr_result.json...")
    eig_csv_path = eigen_sim._output_data_dir + "/eig.csv"
    if not os.path.exists(eig_csv_path):
        print(f"\n!!! FAILED: run() completed without raising, but {eig_csv_path} "
              "was never written (MPI job likely killed mid-run - check for "
              "OOM/signal 9 above).")
        _fail("failed", {"reason": "missing eig.csv, likely killed MPI job"})
        sys.exit(1)

    try:
        interface_epr = eigen_sim.retrieve_interface_EPR_data()
    except Exception:
        print("      WARNING: retrieve_interface_EPR_data() failed - continuing without it.")
        traceback.print_exc()
        interface_epr = None

    port_epr = eigen_sim.retrieve_mode_port_EPR()
    mesh_sizes = eigen_sim.retrieve_simulation_sizes()

    mat_mode_port = port_epr['mat_mode_port']   # (num_modes, num_inductive_ports)
    freqs_hz = port_epr['eigenfrequencies']
    loaded_Q = port_epr['loaded_Q']

    if mat_mode_port is None or mat_mode_port.shape[1] < 1:
        print(f"\n!!! FAILED: port-EPR.csv yielded no inductive-port columns "
              f"(shape {None if mat_mode_port is None else mat_mode_port.shape}).")
        _fail("no_port_epr", {"mat_mode_port_shape": None if mat_mode_port is None
                              else list(mat_mode_port.shape)})
        sys.exit(1)

    n_cols = mat_mode_port.shape[1]
    if n_cols != n_ports:
        print(f"      WARNING: got {n_cols} inductive-port columns, expected "
              f"{n_ports}. Summing over whatever is present.")

    num_modes = len(freqs_hz)
    if num_modes == 0:
        print(f"\n!!! FAILED: zero converged modes. starting_freq="
              f"{args.starting_freq/1e9:.3f}GHz, number_of_freqs={args.number_of_freqs}.")
        _fail("no_modes_converged", {"starting_freq_hz": args.starting_freq,
                                     "number_of_freqs": args.number_of_freqs})
        sys.exit(1)

    import numpy as np
    from build_design import identify_resonator_mode

    p_array_signed = mat_mode_port.sum(axis=1)          # signed total participation
    p_array_abs = np.abs(mat_mode_port).sum(axis=1)     # phonon (~1) vs non-array (<<1)
    per_loop_available = (n_cols == n_ports)
    if per_loop_available:
        p_per_loop = (mat_mode_port[:, 0::2] + mat_mode_port[:, 1::2])  # (num_modes, n_squids)

    # The array's own "phonon" modes store ~all inductive energy in the
    # junctions (sum|p_j| ~ 1); the readout resonator (field on the CPW) has
    # sum|p_j| << 1e-3. Pick the resonator by that, then label the rest
    # phonons by ascending frequency. (min_ghz omitted: sum|p_j| is enough.)
    try:
        res_idx, res_info = identify_resonator_mode(mat_mode_port, freqs_hz)
    except RuntimeError as e:
        print(f"      WARNING: {e}  -> flagging no mode as resonator.")
        res_idx, res_info = None, {"error": str(e)}

    phonon_rank = {}
    for rank, i in enumerate(sorted((i for i in range(num_modes) if i != res_idx),
                                    key=lambda i: freqs_hz[i])):
        phonon_rank[i] = rank + 1

    modes = []
    for i in range(num_modes):
        m = {
            "m": i + 1,
            "f_ghz": float(freqs_hz[i]) / 1e9,
            "Q": float(loaded_Q[i]),
            "p_array_signed": float(p_array_signed[i]),
            "p_array_abs_sum": float(p_array_abs[i]),
            "assignment": "resonator" if i == res_idx else f"phonon_{phonon_rank[i]}",
        }
        if per_loop_available:
            m["p_array_per_loop"] = [float(x) for x in p_per_loop[i]]
        modes.append(m)

    resonator_mode = None
    if res_idx is not None:
        resonator_mode = {
            "m": res_idx + 1,
            "f_ghz": float(freqs_hz[res_idx]) / 1e9,
            "Q": float(loaded_Q[res_idx]),
            "p_array_signed": float(p_array_signed[res_idx]),
            "p_array_abs_sum": float(p_array_abs[res_idx]),
            "p_per_junction": [float(x) for x in mat_mode_port[res_idx]],
        }
        if per_loop_available:
            resonator_mode["p_per_loop"] = [float(x) for x in p_per_loop[res_idx]]

    print("      Modes found:")
    for mode in modes:
        print(f"        m={mode['m']}: f={mode['f_ghz']:.4f}GHz, Q={mode['Q']:.3e}, "
              f"sum|p_j|={mode['p_array_abs_sum']:.3e}, {mode['assignment']}")

    result = {
        "status": "success",
        "n_squids": n_squids,
        "l_j_nh": L_j_nH,
        "cj_j_ff": cj_j_fF,
        "res_length_mm": res_length_mm,
        "geometry_um": geom,
        "n_inductive_ports": n_cols,
        "resonator_id": res_info,
        "resonator_mode": resonator_mode,
        "modes": modes,
        "mesh": {
            "elements": mesh_sizes.get("MeshElements"),
            "resonator_min_size_um": resonator_bg_min_um,
            "array_min_size_um": args.array_bg_min_size_um,
        },
        "interface_epr": _json_safe(interface_epr),
        "wall_time_s": time.time() - t0,
    }
    out_path = f"{args.outdir}/epr_result.json"
    json.dump(result, open(out_path, "w"), indent=2)

    print(f"\n=== STAGE A COMPLETE in {result['wall_time_s']:.1f}s ===")
    print(f"Result written to {out_path}")


if __name__ == "__main__":
    main()
