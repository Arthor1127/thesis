"""
Local, mesh-ONLY test: builds the design and runs SQDMetal's GMSH mesh
generation via prepare_simulation(), then reports element/node counts
directly from GMSH - WITHOUT launching Palace/MPI at all. This is a
seconds-to-minutes check, not a multi-hour cluster job, so use it to
validate a mesh-parameter change before spending cluster time on it.

Runs entirely on your local machine (not the cluster) - just needs a
local Palace binary path for PALACE_Driven_Simulation to accept (it is
never actually invoked here).

Usage - new split-background strategy (recommended, what you're testing):
    python mesh_test_local.py --ground-pos 5.906 \
        --palace-bin /home/ruiz/repo/spack/opt/spack/linux-zen3/palace-develop-blpbtxe43ne2or3yo7gou6ycmdwdxuzl/bin/palace \
        --resonator-bg-min-size-um 30 --resonator-bg-max-size-um 300

Usage - old single-background strategy (for a direct before/after
element-count comparison against the above, same everything else):
    python mesh_test_local.py --ground-pos 5.906 \
        --palace-bin /home/ruiz/repo/spack/opt/spack/linux-zen3/palace-develop-blpbtxe43ne2or3yo7gou6ycmdwdxuzl/bin/palace \
        --old-style

Usage - export the mesh for visual inspection in GMSH or ParaView, e.g.
to actually look at element quality near the port<->resonator and
strap<->resonator size-ratio transitions (Palace log showed kappa maxing
at 175.7, suggesting some badly-shaped elements are sitting somewhere in
one of those transition zones and may be part of why the solve is slow):
    python mesh_test_local.py --ground-pos 5.906 \
        --palace-bin /home/ruiz/repo/spack/opt/spack/linux-zen3/palace-develop-blpbtxe43ne2or3yo7gou6ycmdwdxuzl/bin/palace \
        --resonator-bg-min-size-um 30 --resonator-bg-max-size-um 300 \
        --export mesh_new.msh
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import build_design, TOTAL_LENGTH_MM


def parse_ground_pos(s):
    if s.lower() == "none":
        return None
    return float(s)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ground-pos", type=parse_ground_pos, default=5.906,
                    help="mm, or 'none' for the ungrounded baseline (default: 5.906, "
                         "matching the cluster test_single_driven.sh case)")
    p.add_argument("--palace-bin", required=True,
                    help="Local Palace binary path. Never actually invoked in this "
                         "script (no .run() call) but PALACE_Driven_Simulation wants "
                         "a real-looking path at construction time.")
    p.add_argument("--total-length-mm", type=float, default=None)
    p.add_argument("--bg-min-size-um", type=float, default=15.0,
                    help="Background mesh min_size (um) for TL/LP1/LP2 (ports). "
                         "In --old-style mode, this also applies to resonator1.")
    p.add_argument("--bg-max-size-um", type=float, default=150.0)
    p.add_argument("--resonator-bg-min-size-um", type=float, default=None,
                    help="resonator1's own background min_size (um). Defaults to "
                         "--bg-min-size-um if omitted. Ignored in --old-style mode.")
    p.add_argument("--resonator-bg-max-size-um", type=float, default=None)
    p.add_argument("--strap-margin-um", type=float, default=3.0)
    p.add_argument("--resonator-min-size-um", type=float, default=2.0,
                    help="Fine mesh size (um) inside the strap-local box only.")
    p.add_argument("--old-style", action="store_true",
                    help="Use the ORIGINAL single fine_mesh_components() call over "
                         "['TL','resonator1','LP1','LP2'] together at --bg-min/max-"
                         "size-um, instead of the new split (ports vs resonator1) "
                         "strategy. Use this to reproduce the old element count for "
                         "direct comparison against a --resonator-bg-* run.")
    p.add_argument("--export", default=None,
                    help="Path to write the generated mesh for visual inspection, "
                         "e.g. mesh_new.msh (open in GMSH) or mesh_new.vtk (open in "
                         "ParaView). Format is inferred from the extension - GMSH "
                         "supports .msh/.vtk/.vtu/.stl and others directly via "
                         "gmsh.write(). Useful for actually LOOKING at element "
                         "quality/kappa near the port<->resonator and strap<->"
                         "resonator size-ratio transitions, since the Palace log's "
                         "kappa=175.7 suggests some badly-shaped elements are sitting "
                         "somewhere in exactly one of those transition zones.")
    args = p.parse_args()

    resonator_bg_min_um = (args.bg_min_size_um if args.old_style else
                            (args.resonator_bg_min_size_um
                             if args.resonator_bg_min_size_um is not None
                             else args.bg_min_size_um))
    resonator_bg_max_um = (args.bg_max_size_um if args.old_style else
                            (args.resonator_bg_max_size_um
                             if args.resonator_bg_max_size_um is not None
                             else args.bg_max_size_um))

    ground_pos_label = "none" if args.ground_pos is None else f"{args.ground_pos}mm"
    total_length_mm = args.total_length_mm if args.total_length_mm is not None else TOTAL_LENGTH_MM
    strategy = "OLD (single background, resonator1 == ports)" if args.old_style else "NEW (split background)"

    print(f"=== Mesh-only test: strategy={strategy}, ground_pos={ground_pos_label} ===")
    print(f"    ports (TL/LP1/LP2): {args.bg_min_size_um}/{args.bg_max_size_um}um")
    print(f"    resonator1:         {resonator_bg_min_um}/{resonator_bg_max_um}um")
    print(f"    strap-local box:    {args.resonator_min_size_um}um "
          f"(margin={args.strap_margin_um}um)")

    t0 = time.time()
    print("\n[1/3] Building design...")
    design = build_design(ground_pos_mm=args.ground_pos, total_length_mm=total_length_mm)
    design.rebuild()
    print(f"      OK ({time.time()-t0:.1f}s)")

    print("[2/3] Setting up driven simulation object (Palace NOT invoked)...")
    from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation
    import matplotlib
    matplotlib.use("Agg", force=True)

    driven_options = {
        "dielectric_material": "silicon",
        "solns_to_save": 1,
        "solver_order": 2,
        "solver_tol": 1.0e-8,
        "solver_maxits": 500,
        "fillet_resolution": 12,
        "palace_dir": args.palace_bin,
        "num_cpus": 1,
    }

    outdir = os.path.expanduser("~/mesh_test_scratch")
    os.makedirs(outdir, exist_ok=True)

    driven_sim = PALACE_Driven_Simulation(
        name=f"meshtest_{ground_pos_label}",
        metal_design=design,
        sim_parent_directory=outdir,
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

    if args.old_style:
        driven_sim.fine_mesh_components(
            ["TL", "resonator1", "LP1", "LP2"],
            min_size=args.bg_min_size_um * 1e-6, max_size=args.bg_max_size_um * 1e-6,
            taper_dist_min=10e-6, metals_only=False,
        )
    else:
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
            print("WARNING: --ground-pos given but no ground_strap_* geometry found.")
        else:
            strap_geom = strap_rows.iloc[0]["geometry"]
            minx, miny, maxx, maxy = strap_geom.bounds
            safe_margin_mm = args.strap_margin_um / 1000
            x1, y1 = (minx - safe_margin_mm) * 1e-3, (miny - safe_margin_mm) * 1e-3
            x2, y2 = (maxx + safe_margin_mm) * 1e-3, (maxy + safe_margin_mm) * 1e-3
            taper_target_um = args.bg_min_size_um if args.old_style else resonator_bg_min_um
            driven_sim.fine_mesh_in_rectangle(
                x1, y1, x2, y2,
                min_size=args.resonator_min_size_um * 1e-6, max_size=taper_target_um * 1e-6,
                taper_dist_min=2e-6, taper_dist_max=6e-6,
            )

    driven_sim.set_freq_values(freq_start=5e9, freq_end=5e9, freq_step=1e9)

    print("[3/3] Running prepare_simulation() - this triggers GMSH mesh generation "
          "only, no Palace/MPI launch...")
    t1 = time.time()
    driven_sim.prepare_simulation()
    print(f"      Mesh generated in {time.time()-t1:.1f}s")

    import gmsh
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements()
    total_elements = sum(len(t) for t in elem_tags)

    if args.export:
        export_path = os.path.abspath(os.path.expanduser(args.export))
        gmsh.write(export_path)
        print(f"\n    Mesh written to: {export_path}")
        print(f"    Open directly in GMSH (File > Open), or in ParaView if .vtk/.vtu.")

    print(f"\n=== RESULT ({strategy}) ===")
    print(f"    Nodes:    {len(node_tags):,}")
    print(f"    Elements: {total_elements:,}")
    print(f"    Total wall time: {time.time()-t0:.1f}s")
    print("\nCompare this Elements count against a run with the opposite "
          "--old-style setting (same --ground-pos) to see the actual effect "
          "of the mesh split before spending cluster time on it.")


if __name__ == "__main__":
    main()
