"""
One step of a fine-mesh-size sweep: build the design, apply fine meshing
to ONE region (readout / pad / squid_arm / junction) at a swept min_size
while the other three regions stay at a coarse baseline, then call
prepare_simulation() (GMSH meshing + Palace config write, NO Palace
solve). Meant to be run under `timeout 120 python -u mesh_sweep_step.py
...` so a step that hangs in GMSH's background-field computation gets
killed cleanly by the wrapper rather than this script.

On success, prints "MESH_OK elements=<n> wall_s=<t>" and writes a JSON
result file. On any exception, prints "MESH_FAILED reason=<...>".
(A timeout is NOT caught here - it shows up as no result file + the
`timeout` wrapper's own nonzero exit, which the caller checks for.)
"""
import argparse
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import (build_design, L_TRANSMON_NH, L_PROBE_NH,
                           CJ_TRANSMON_FF, CJ_PROBE_FF, QUBIT_COMPONENT_NAME,
                           JJ_INDEX_SINGLE, JJ_INDEX_SQUID_A)

BASELINE_UM = 25.0  # coarse-test value, used for the three regions NOT being swept

# Squid-arm bbox computed from squid_options geometry (theta_2=180,
# delta_angle=7.5, pad_radius=90um, l1=l2=300um, w1=6um, g1=4um) - see
# 2026-08-24 chat derivation. Fixed regardless of which region is swept.
SQUID_ARM_BOX_M = (-710e-6, -25e-6, -80e-6, 25e-6)
PAD_RADIUS_M = 90e-6
PAD_MARGIN_M = 20e-6


def max_size_for(min_size_um, cap_um):
    return min(max(min_size_um * 4, 20.0), cap_um) * 1e-6


def count_msh_elements(msh_path):
    with open(msh_path) as f:
        for line in f:
            if line.strip() == "$Elements":
                return int(next(f).strip())
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", choices=["readout", "pad", "squid_arm", "junction"])
    ap.add_argument("--min-size-um", type=float)
    ap.add_argument("--all-min-size-um", type=float, default=None,
                     help="If given, sets ALL four regions to this value "
                          "(instead of one region swept vs. the others at baseline). "
                          "--region/--min-size-um are ignored in that case.")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    if args.all_min_size_um is not None:
        sizes_um = {k: args.all_min_size_um for k in
                    ("readout", "pad", "squid_arm", "junction")}
        args.region = f"all@{args.all_min_size_um}um"
        args.min_size_um = args.all_min_size_um
    else:
        if args.region is None or args.min_size_um is None:
            ap.error("--region and --min-size-um are required unless --all-min-size-um is given")
        sizes_um = {"readout": BASELINE_UM, "pad": BASELINE_UM,
                    "squid_arm": BASELINE_UM, "junction": BASELINE_UM}
        sizes_um[args.region] = args.min_size_um

    os.makedirs(args.outdir, exist_ok=True)
    result_path = os.path.join(args.outdir, "result.json")

    t0 = time.time()
    try:
        design = build_design()
        design.rebuild()

        from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
        from SQDMetal.Utilities.Materials import MaterialInterface
        import matplotlib
        matplotlib.use("Agg", force=True)

        user_defined_options = {
            "mesh_refinement": 0,
            "dielectric_material": "silicon",
            "starting_freq": 7.5e9,
            "number_of_freqs": 6,
            "solns_to_save": 6,
            "solver_order": 2,
            "solver_tol": 1.0e-8,
            "solver_maxits": 200,
            "fillet_resolution": 12,
            "palace_dir": os.environ.get("PALACE_BIN", "/nonexistent"),
            "num_cpus": 1,
        }
        eigen_sim = PALACE_Eigenmode_Simulation(
            name="sweep_step",
            metal_design=design,
            sim_parent_directory=args.outdir,
            mode="simPC",
            meshing="GMSH",
            user_options=user_defined_options,
            create_files=True,
        )
        eigen_sim.add_metallic(1)
        eigen_sim.add_ground_plane()
        eigen_sim.create_port_JosephsonJunction(
            QUBIT_COMPONENT_NAME, junction_index=JJ_INDEX_SINGLE,
            L_J=L_TRANSMON_NH * 1e-9, C_J=CJ_TRANSMON_FF * 1e-15)
        eigen_sim.create_port_JosephsonJunction(
            QUBIT_COMPONENT_NAME, junction_index=JJ_INDEX_SQUID_A,
            L_J=L_PROBE_NH * 1e-9, C_J=CJ_PROBE_FF * 1e-15)

        eigen_sim.fine_mesh_components(
            ["readout_1"],
            min_size=sizes_um["readout"] * 1e-6, max_size=max_size_for(sizes_um["readout"], 150),
            taper_dist_min=10e-6, metals_only=False,
        )
        eigen_sim.fine_mesh_in_rectangle(
            -(PAD_RADIUS_M + PAD_MARGIN_M), -(PAD_RADIUS_M + PAD_MARGIN_M),
            (PAD_RADIUS_M + PAD_MARGIN_M), (PAD_RADIUS_M + PAD_MARGIN_M),
            min_size=sizes_um["pad"] * 1e-6, max_size=max_size_for(sizes_um["pad"], 150),
            taper_dist_min=5e-6, taper_dist_max=20e-6,
        )
        eigen_sim.fine_mesh_in_rectangle(
            *SQUID_ARM_BOX_M,
            min_size=sizes_um["squid_arm"] * 1e-6, max_size=max_size_for(sizes_um["squid_arm"], 150),
            taper_dist_min=5e-6, taper_dist_max=20e-6,
        )
        for port in eigen_sim._ports[-2:]:
            coords = port['portCoords']
            xs = [c[0] for c in coords]
            ys = [c[1] for c in coords]
            eigen_sim.fine_mesh_in_rectangle(
                min(xs), min(ys), max(xs), max(ys),
                min_size=sizes_um["junction"] * 1e-6, max_size=max_size_for(sizes_um["junction"], 50),
                taper_dist_min=5e-6, taper_dist_max=20e-6,
            )
        eigen_sim.setup_EPR_interfaces(
            metal_air=MaterialInterface("Aluminium-Vacuum"),
            substrate_air=MaterialInterface("Silicon-Vacuum"),
            substrate_metal=MaterialInterface("Silicon-Aluminium"),
        )
        eigen_sim.prepare_simulation()
        wall_s = time.time() - t0

        msh_path = os.path.join(args.outdir, "sweep_step", "sweep_step.msh")
        elements = count_msh_elements(msh_path) if os.path.exists(msh_path) else None

        result = {"status": "ok", "region": args.region, "min_size_um": args.min_size_um,
                  "sizes_um": sizes_um, "wall_s": wall_s, "elements": elements}
        json.dump(result, open(result_path, "w"), indent=2)
        print(f"MESH_OK elements={elements} wall_s={wall_s:.1f}", flush=True)
    except Exception as e:
        wall_s = time.time() - t0
        result = {"status": "failed", "region": args.region, "min_size_um": args.min_size_um,
                  "sizes_um": sizes_um, "wall_s": wall_s, "reason": str(e)}
        json.dump(result, open(result_path, "w"), indent=2)
        print(f"MESH_FAILED reason={e}", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
