"""Palace eigenmode ("normal modes") solve for the bare Josephson-junction array.

Builds the design from ``build_design.py``, puts one lumped ``L_j``/``C_j``
inductive port across every junction gap, and runs a Palace eigenmode solve.
Writes ``modes.json`` with the eigenfrequencies, loaded Q, and the per-junction
participation vector of each mode (handy for seeing which standing-wave pattern
each normal mode is).

NO EPR / dielectric-loss analysis is done (``setup_EPR_interfaces`` is not
called, so no surface-Q.csv) - the user asked for the eigenmodes only. The
per-junction participation still comes for free from Palace's ``port-EPR.csv``.

Adapted from ``../squid_array_coupling/build_and_eigenmode.py`` (drop the
resonator, the coupling-zone boxes, and the interface EPR).
"""

import argparse
import json
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import (build_design, N_JUNCTIONS, L_J_NH, CJ_J_FF, JJ_WIDTH_UM,
                          PAD_LENGTH_UM, PAD_WIDTH_UM, PAD_GAP_UM, GROUND_GAP_UM,
                          ARRAY_COMPONENT_NAME)


def _build_parser():
    p = argparse.ArgumentParser(
        description="Build the bare JJ array and run a Palace eigenmode solve. "
                    "Any knob defaulting to None falls back to the constant of "
                    "the same name in build_design.py.")

    sweep = p.add_argument_group("geometry / inductance")
    sweep.add_argument("--n-junctions", type=int, default=None,
                       help=f"Junctions in the chain. Default {N_JUNCTIONS}.")
    sweep.add_argument("--L-j-nH", type=float, default=None,
                       help=f"Inductance of every junction (nH). Default {L_J_NH}.")
    sweep.add_argument("--Cj-j-fF", type=float, default=None,
                       help=f"Capacitance of every junction (fF). Default {CJ_J_FF}.")
    sweep.add_argument("--jj-width-um", type=float, default=None,
                       help=f"Junction bridge transverse width (um). Default {JJ_WIDTH_UM}.")
    sweep.add_argument("--pad-length-um", type=float, default=None,
                       help=f"Island length along the chain (um). Default {PAD_LENGTH_UM}.")
    sweep.add_argument("--pad-width-um", type=float, default=None,
                       help=f"Island transverse width (um). Default {PAD_WIDTH_UM}.")
    sweep.add_argument("--pad-gap-um", type=float, default=None,
                       help=f"Vacuum break each junction bridges (um). Default {PAD_GAP_UM}.")
    sweep.add_argument("--ground-gap-um", type=float, default=None,
                       help=f"Ground-plane pocket width around the chain (um). Default {GROUND_GAP_UM}.")

    solver = p.add_argument_group("eigenmode solver")
    solver.add_argument("--number-of-freqs", type=int, default=8,
                        help="Number of eigenmodes to solve for.")
    solver.add_argument("--starting-freq", type=float, default=1e9,
                        help="Eigen-search start frequency (Hz). The array's "
                             "fundamental plasma mode is the lowest; start below it.")
    solver.add_argument("--solver-order", type=int, default=1)
    solver.add_argument("--solver-tol", type=float, default=1.0e-5)
    solver.add_argument("--solver-maxits", type=int, default=300)
    solver.add_argument("--fillet-resolution", type=int, default=12)

    mesh = p.add_argument_group("mesh")
    mesh.add_argument("--bg-min-size-um", type=float, default=8.0,
                      help="General background min element size (um).")
    mesh.add_argument("--bg-max-size-um", type=float, default=120.0)
    mesh.add_argument("--array-bg-min-size-um", type=float, default=2.0,
                      help="Array fine-mesh min size (um). Finest features are the "
                           "pad_gap (~0.15um, the junction port) and the island "
                           "edges; 2um is the practical GMSH floor seen in the "
                           "SQUID-array runs.")
    mesh.add_argument("--array-bg-max-size-um", type=float, default=120.0)
    mesh.add_argument("--array-taper-dist-min-um", type=float, default=10.0,
                      help="Distance from the array metal held at min size (um).")
    mesh.add_argument("--array-taper-dist-max-um", type=float, default=200.0)

    chip = p.add_argument_group("chip / infra")
    chip.add_argument("--chip-margin-um", type=float, default=None)
    chip.add_argument("--chip-size-x-mm", type=float, default=None)
    chip.add_argument("--chip-size-y-mm", type=float, default=None)
    chip.add_argument("--chip-center-x-mm", type=float, default=None)
    chip.add_argument("--chip-center-y-mm", type=float, default=None)
    chip.add_argument("--outdir", required=True)
    chip.add_argument("--palace-bin", required=True)
    chip.add_argument("--num-cpus", type=int, default=32)

    return p


def main():
    args = _build_parser().parse_args()

    def _d(v, default):
        return default if v is None else v

    n_junctions = _d(args.n_junctions, N_JUNCTIONS)
    L_j_nH = _d(args.L_j_nH, L_J_NH)
    cj_j_fF = _d(args.Cj_j_fF, CJ_J_FF)

    geom = dict(
        jj_width_um=_d(args.jj_width_um, JJ_WIDTH_UM),
        pad_length_um=_d(args.pad_length_um, PAD_LENGTH_UM),
        pad_width_um=_d(args.pad_width_um, PAD_WIDTH_UM),
        pad_gap_um=_d(args.pad_gap_um, PAD_GAP_UM),
        ground_gap_um=_d(args.ground_gap_um, GROUND_GAP_UM),
    )

    os.makedirs(args.outdir, exist_ok=True)

    def _fail(status, extra=None):
        payload = {"status": status, "n_junctions": n_junctions,
                   "l_j_nh": L_j_nH, "cj_j_ff": cj_j_fF}
        if extra:
            payload.update(extra)
        json.dump(payload, open(f"{args.outdir}/modes.json", "w"), indent=2)

    t0 = time.time()
    print(f"[1/4] Building design (n_junctions={n_junctions}, L_j={L_j_nH}nH, "
          f"C_j={cj_j_fF}fF, ground_gap={geom['ground_gap_um']}um)...")
    build_kwargs = dict(geom)
    if args.chip_margin_um is not None:
        build_kwargs["chip_margin_um"] = args.chip_margin_um
    design = build_design(n_junctions=n_junctions, L_j_nH=L_j_nH, cj_j_fF=cj_j_fF,
                          chip_size_x_mm=args.chip_size_x_mm,
                          chip_size_y_mm=args.chip_size_y_mm,
                          chip_center_x_mm=args.chip_center_x_mm,
                          chip_center_y_mm=args.chip_center_y_mm,
                          **build_kwargs)
    design.rebuild()

    n_jj = len(design.qgeometry.tables['junction'].loc[
        design.qgeometry.tables['junction']['component']
        == design.components[ARRAY_COMPONENT_NAME].id])
    print(f"      OK ({time.time()-t0:.1f}s). Components: {list(design.components.keys())}; "
          f"{n_jj} junction rows.")
    if n_jj != n_junctions:
        print(f"\n!!! FAILED: {n_jj} junction rows, expected {n_junctions}.")
        _fail("failed", {"reason": "junction count mismatch", "n_junction_rows": n_jj})
        sys.exit(1)

    print("[2/4] Importing SQDMetal PALACE eigenmode module...")
    from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
    print("      OK")

    import matplotlib
    matplotlib.use("Agg", force=True)

    print(f"[3/4] Setting up + running eigenmode solve "
          f"({n_junctions} lumped inductive ports on the array)...")
    t1 = time.time()
    user_defined_options = {
        "mesh_refinement": 0,
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

    run_label = f"njj_{n_junctions}_Lj_{L_j_nH}nH"

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

        # One lumped inductive port per junction gap. JosephsonArray emits the
        # junctions in order jj_0 .. jj_{n-1}; every gap is a real vacuum break
        # and MUST be ported or the chain is open there.
        for k in range(n_junctions):
            eigen_sim.create_port_JosephsonJunction(
                ARRAY_COMPONENT_NAME, junction_index=k,
                L_J=L_j_nH * 1e-9, C_J=cj_j_fF * 1e-15)

        eigen_sim.fine_mesh_components(
            [ARRAY_COMPONENT_NAME],
            min_size=args.array_bg_min_size_um * 1e-6,
            max_size=args.array_bg_max_size_um * 1e-6,
            taper_dist_min=args.array_taper_dist_min_um * 1e-6,
            taper_dist_max=args.array_taper_dist_max_um * 1e-6,
            metals_only=False,
        )

        # No setup_EPR_interfaces(): eigenmodes only, no dielectric-loss EPR.

        eigen_sim.prepare_simulation()
        eigen_sim.run()
    except Exception:
        print("\n!!! FAILED during eigenmode simulation !!!")
        traceback.print_exc()
        _fail("failed")
        sys.exit(1)

    print(f"      OK ({time.time()-t1:.1f}s)")

    print("[4/4] Retrieving eigenmodes and writing modes.json...")
    eig_csv_path = eigen_sim._output_data_dir + "/eig.csv"
    if not os.path.exists(eig_csv_path):
        print(f"\n!!! FAILED: run() returned but {eig_csv_path} was never written "
              "(MPI job likely killed mid-run - check for OOM/signal 9 above).")
        _fail("failed", {"reason": "missing eig.csv, likely killed MPI job"})
        sys.exit(1)

    import numpy as np

    port_epr = eigen_sim.retrieve_mode_port_EPR()
    mesh_sizes = eigen_sim.retrieve_simulation_sizes()

    mat_mode_port = port_epr['mat_mode_port']   # (num_modes, num_inductive_ports)
    freqs_hz = port_epr['eigenfrequencies']
    loaded_Q = port_epr['loaded_Q']

    num_modes = len(freqs_hz)
    if num_modes == 0:
        print(f"\n!!! FAILED: zero converged modes. starting_freq="
              f"{args.starting_freq/1e9:.3f}GHz, number_of_freqs={args.number_of_freqs}.")
        _fail("no_modes_converged", {"starting_freq_hz": args.starting_freq,
                                     "number_of_freqs": args.number_of_freqs})
        sys.exit(1)

    have_ports = mat_mode_port is not None and mat_mode_port.shape[1] >= 1

    modes = []
    for i in range(num_modes):
        m = {
            "m": i + 1,
            "f_ghz": float(freqs_hz[i]) / 1e9,
            "Q": float(loaded_Q[i]),
        }
        if have_ports:
            p_j = mat_mode_port[i]
            m["p_array_abs_sum"] = float(np.abs(p_j).sum())
            m["p_per_junction"] = [float(x) for x in p_j]
        modes.append(m)

    print("      Modes found:")
    for mode in modes:
        extra = (f", sum|p_j|={mode['p_array_abs_sum']:.3e}"
                 if "p_array_abs_sum" in mode else "")
        print(f"        m={mode['m']}: f={mode['f_ghz']:.5f}GHz, Q={mode['Q']:.3e}{extra}")

    result = {
        "status": "success",
        "n_junctions": n_junctions,
        "l_j_nh": L_j_nH,
        "cj_j_ff": cj_j_fF,
        "geometry_um": geom,
        "n_inductive_ports": (mat_mode_port.shape[1] if have_ports else 0),
        "chip_mm": {
            "size_x": design.chips.main.size.size_x,
            "size_y": design.chips.main.size.size_y,
            "size_z": design.chips.main.size.size_z,
        },
        "mesh": {
            "elements": mesh_sizes.get("MeshElements"),
            "array_min_size_um": args.array_bg_min_size_um,
        },
        "modes": modes,
        "wall_time_s": time.time() - t0,
    }
    out_path = f"{args.outdir}/modes.json"
    json.dump(result, open(out_path, "w"), indent=2)

    print(f"\n=== EIGENMODE SOLVE COMPLETE in {result['wall_time_s']:.1f}s ===")
    print(f"Result written to {out_path}")


if __name__ == "__main__":
    main()
