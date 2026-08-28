import argparse
import json
import os
import sys
import time
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import (build_design, L_TRANSMON_NH, L_PROBE_NH, RES_LENGTH_MM,
                           CJ_TRANSMON_FF, CJ_PROBE_FF, QUBIT_COMPONENT_NAME,
                           JJ_INDEX_SINGLE, JJ_INDEX_SQUID_A)


def _json_safe(obj):
    """Recursively convert numpy arrays/scalars into plain Python types so
    json.dump doesn't choke on them.

    numpy complex scalars have .tolist(), but it converts them to a plain
    Python `complex` - which is STILL not JSON serializable, so that result
    must be re-checked (2026-08-25: this crashed mid-write on interface_epr
    data, corrupting the output file for a run whose actual solve had
    already succeeded)."""
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
    p.add_argument("--L-transmon-nH", type=float, default=None,
                    help="Transmon junction (to ground) inductance L_T (nH). "
                        "Defaults to build_design.L_TRANSMON_NH.")
    p.add_argument("--L-probe-nH", type=float, default=None,
                    help="SQUID branch-a probe inductance L_probe (nH). "
                        "Defaults to build_design.L_PROBE_NH.")
    p.add_argument("--res-length-mm", type=float, default=None,
                    help="Readout resonator total length (mm). Defaults to build_design.RES_LENGTH_MM.")
    p.add_argument("--outdir", required=True)
    p.add_argument("--palace-bin", required=True)
    p.add_argument("--num-cpus", type=int, default=16)
    p.add_argument("--number-of-freqs", type=int, default=6,
                    help="Number of eigenmodes to solve. EPR_C_EXTRACTION.md Sec 4.3 "
                        "recommends at least 6 (resonator fundamental + harmonics + transmon).")
    p.add_argument("--starting-freq", type=float, default=5e9)
    p.add_argument("--bg-min-size-um", type=float, default=6.0,
                    help="General background min element size (um); also the fallback "
                        "for --resonator-bg-min-size-um. 6um resolves trace_gap=12.25um.")
    p.add_argument("--bg-max-size-um", type=float, default=150.0)
    p.add_argument("--amr-max-its", type=int, default=0)
    p.add_argument("--amr-tol", type=float, default=1e-2)
    p.add_argument("--chip-size-x-mm", type=float, default=None)
    p.add_argument("--chip-size-y-mm", type=float, default=None)
    p.add_argument("--chip-center-x-mm", type=float, default=None)
    p.add_argument("--chip-center-y-mm", type=float, default=None)
    p.add_argument("--resonator-bg-min-size-um", type=float, default=None,
                    help="readout_1's own background min_size (um). Should resolve "
                        "trace_gap=12.25um (design.variables['cpw_gap']) - roughly "
                        "half the gap, i.e. ~6um, is a reasonable floor. Defaults to "
                        "--bg-min-size-um if not given.")
    p.add_argument("--resonator-bg-max-size-um", type=float, default=None,
                    help="readout_1's own background max_size (um). "
                        "Defaults to --bg-max-size-um if not given.")
    p.add_argument("--qubit-bg-min-size-um", type=float, default=2.0,
                    help="Q1_SQUID's own background min_size (um), separate from the "
                        "resonator's since the qubit has finer features. Smallest "
                        "non-junction-port gap on the qubit is squid_options.g1=4um "
                        "(SQUID branch-to-ground clearance); ~half that, 2um, "
                        "resolves it the same way --resonator-bg-min-size-um resolves "
                        "trace_gap. (pad_gap=10um is coarser and already covered by "
                        "this.) Junction gaps themselves (jj_sim_gap=3um) are handled "
                        "separately by --junction-mesh-min-um.")
    p.add_argument("--qubit-bg-max-size-um", type=float, default=150.0,
                    help="Q1_SQUID's own background max_size (um). Raised from an "
                        "original 20.0um default (2026-08-25): 20um forced a large "
                        "fraction of the whole domain to stay at that size instead of "
                        "relaxing further, a major contributor to an 11.5M-element "
                        "mesh blowup at otherwise-identical settings.")
    p.add_argument("--junction-mesh-min-um", type=float, default=6.0,
                    help="Fine-mesh min element size local to each junction port (um). "
                        "Raised from an original 1.5um default per the 2026-08-25 "
                        "per-region sweep: EVERY region tested (readout, pad, squid_arm, "
                        "junction) independently passed at 6um (~80-96s) and timed out "
                        "(>120s) at 5um and below - a shared GMSH background-field-"
                        "computation floor, not something specific to one region's size "
                        "or shape. Do not go below 6um without re-validating against that "
                        "floor first.")
    p.add_argument("--junction-mesh-max-um", type=float, default=20.0,
                    help="Fine-mesh max element size local to each junction port (um).")
    args = p.parse_args()

    L_transmon_nH = args.L_transmon_nH if args.L_transmon_nH is not None else L_TRANSMON_NH
    L_probe_nH = args.L_probe_nH if args.L_probe_nH is not None else L_PROBE_NH
    res_length_mm = args.res_length_mm if args.res_length_mm is not None else RES_LENGTH_MM

    t0 = time.time()
    print(f"[1/4] Building design (L_transmon={L_transmon_nH}nH, L_probe={L_probe_nH}nH, "
          f"res_length={res_length_mm}mm)...")
    design = build_design(L_transmon_nH=L_transmon_nH, L_probe_nH=L_probe_nH,
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

    print("[3/4] Setting up + running eigenmode simulation (production settings)...")
    t1 = time.time()
    user_defined_options = {
        "mesh_refinement": 0,  # UniformLevels - NOT AMR; see --amr-max-its
        "dielectric_material": "silicon",
        "starting_freq": args.starting_freq,
        "number_of_freqs": args.number_of_freqs,
        "solns_to_save": args.number_of_freqs,
        # order=1/tol=1e-5 (not the original order=2/tol=1e-8) per the
        # 2026-08-25 findings: order=2 roughly 6-8x's the DOF count for the
        # same mesh and was the direct cause of an earlier 38.6M-unknown
        # OOM (the "6um-everywhere" run); order=1 has converged cleanly
        # (1 FGMRES iteration per solve) at every mesh size tested since.
        "solver_order": 1,
        "solver_tol": 1.0e-5,
        "solver_maxits": 300,
        "fillet_resolution": 12,
        "palace_dir": args.palace_bin,
        "num_cpus": args.num_cpus,
    }

    run_label = f"Ltransmon_{L_transmon_nH}_Lprobe_{L_probe_nH}_reslen_{res_length_mm}mm"

    def _fail(status, extra=None):
        payload = {"status": status, "L_transmon_nH": L_transmon_nH,
                   "L_probe_nH": L_probe_nH, "res_length_mm": res_length_mm}
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

        # Exactly two lumped inductive ports, per EPR_C_EXTRACTION.md Sec 4.1:
        # the transmon junction (L_T) and the SQUID's probe arm a (L_probe).
        # Branch b is bridged with continuous metal in build_design.py
        # (jj_b_options.export_mask=True) and gets no port at all - it is not
        # a real junction for this measurement. port-EPR.csv will therefore
        # have exactly 2 participation columns, in this creation order:
        # column 0 = transmon, column 1 = probe. SQDMetal has no name-based
        # port lookup (create_port_JosephsonJunction only stores an
        # auto-generated port_name, not the source component/junction_index),
        # so this ordering must be tracked here rather than read back later.
        eigen_sim.create_port_JosephsonJunction(
            QUBIT_COMPONENT_NAME, junction_index=JJ_INDEX_SINGLE,
            L_J=L_transmon_nH * 1e-9, C_J=CJ_TRANSMON_FF * 1e-15)
        eigen_sim.create_port_JosephsonJunction(
            QUBIT_COMPONENT_NAME, junction_index=JJ_INDEX_SQUID_A,
            L_J=L_probe_nH * 1e-9, C_J=CJ_PROBE_FF * 1e-15)

        resonator_bg_min_um = (args.resonator_bg_min_size_um
                                if args.resonator_bg_min_size_um is not None
                                else args.bg_min_size_um)
        resonator_bg_max_um = (args.resonator_bg_max_size_um
                                if args.resonator_bg_max_size_um is not None
                                else args.bg_max_size_um)
        # taper_dist_min widened per the 2026-08-25 AMR run: the Indicator
        # field showed refinement concentrating broadly around the readout
        # resonator (mode 1's field lives entirely on its trace) while the
        # qubit pad/SQUID arms as a whole only accounted for ~1-4% of
        # indicator mass - so the readout's fine zone extends further and
        # the qubit's stays narrow.
        eigen_sim.fine_mesh_components(
            ["readout_1"],
            min_size=resonator_bg_min_um * 1e-6, max_size=resonator_bg_max_um * 1e-6,
            taper_dist_min=40e-6, metals_only=False,
        )
        eigen_sim.fine_mesh_components(
            [QUBIT_COMPONENT_NAME],
            min_size=args.qubit_bg_min_size_um * 1e-6, max_size=args.qubit_bg_max_size_um * 1e-6,
            taper_dist_min=5e-6, metals_only=False,
        )

        # A dedicated squid-loop fine_mesh_in_rectangle was tried here and
        # reverted (2026-08-25): fine_mesh_in_rectangle fine-meshes the
        # WHOLE box regardless of whether there's metal there, unlike
        # fine_mesh_components (which only refines near actual geometry
        # edges) - most of a box wide enough to cover the squid arms is
        # empty substrate/ground plane, so it blew unknowns up from 2.86M
        # to 13.2M for no benefit. A direct mesh-file check of the
        # qubit-wide pass alone already showed ~3.1um median element size
        # right at the squid loop, so it was already adequately resolved
        # without this addition.

        # Junction-local fine_mesh_in_rectangle REMOVED (2026-08-25): an
        # isolation test (identical config with/without this block) showed
        # it alone - not the readout taper widening, not the qubit max_size
        # - was what pushed mesh generation from 52.9s/2.47M elements to
        # a >180s timeout, even though the port bounding box itself is
        # tiny. Root cause not fully understood (some background-field
        # interaction with the other fine-mesh regions), but the qubit-wide
        # pass alone already gives ~3.1um median element size at the SQUID
        # loop (verified directly from a generated mesh file), so junction
        # resolution was not actually depending on this block.
        # --junction-mesh-min-um/--junction-mesh-max-um are now unused.

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
              "was never written.")
        print("This usually means the MPI job itself was killed mid-run "
              "(check for OOM/signal 9 in the log above) without Palace/SQDMetal "
              "propagating that as a Python exception.")
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

    mat_mode_port = port_epr['mat_mode_port']  # shape (num_modes, num_inductive_ports)
    freqs_hz = port_epr['eigenfrequencies']
    loaded_Q = port_epr['loaded_Q']

    if mat_mode_port is None or mat_mode_port.shape[1] < 2:
        print(f"\n!!! FAILED: port-EPR.csv did not yield 2 inductive-port columns "
              f"(got shape {None if mat_mode_port is None else mat_mode_port.shape}). "
              "Check that both create_port_JosephsonJunction calls actually produced "
              "inductive ports (L_J > 0).")
        _fail("no_port_epr", {"mat_mode_port_shape": None if mat_mode_port is None
                               else list(mat_mode_port.shape)})
        sys.exit(1)

    num_modes = len(freqs_hz)
    if num_modes == 0:
        print(f"\n!!! FAILED: zero converged modes. starting_freq="
              f"{args.starting_freq/1e9:.3f}GHz, number_of_freqs={args.number_of_freqs}.")
        _fail("no_modes_converged", {"starting_freq_hz": args.starting_freq,
                                      "number_of_freqs": args.number_of_freqs})
        sys.exit(1)

    p_transmon = mat_mode_port[:, 0]
    p_probe = mat_mode_port[:, 1]

    # Participation-based mode ID (EPR_C_EXTRACTION.md Sec 6.3): sort on
    # participation, not frequency - frequency-based ID fails silently the
    # moment a qubit mode crosses a resonator harmonic.
    transmon_idx = int(p_transmon.argmax())
    freq_order = sorted((i for i in range(num_modes) if i != transmon_idx),
                         key=lambda i: freqs_hz[i])
    resonator_rank = {i: rank + 1 for rank, i in enumerate(freq_order)}

    modes = []
    for i in range(num_modes):
        assignment = "transmon" if i == transmon_idx else f"resonator_mode_{resonator_rank[i]}"
        modes.append({
            "m": i + 1,
            "f_ghz": float(freqs_hz[i]) / 1e9,
            "Q": float(loaded_Q[i]),
            "p": {"probe": float(p_probe[i]), "transmon": float(p_transmon[i])},
            "assignment": assignment,
        })

    print("      Modes found:")
    for mode in modes:
        print(f"        m={mode['m']}: f={mode['f_ghz']:.4f}GHz, Q={mode['Q']:.3e}, "
              f"p_probe={mode['p']['probe']:.3e}, p_transmon={mode['p']['transmon']:.3e}, "
              f"assignment={mode['assignment']}")

    result = {
        "status": "success",
        "l_probe_nh": L_probe_nH,
        "l_transmon_nh": L_transmon_nH,
        "res_length_mm": res_length_mm,
        "modes": modes,
        "mesh": {
            "elements": mesh_sizes.get("MeshElements"),
            "min_size_um": resonator_bg_min_um,
            "qubit_min_size_um": args.qubit_bg_min_size_um,
            "junction_min_size_um": args.junction_mesh_min_um,
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
