"""
Local, mesh-ONLY test for the quarter-wave design: builds the design and
runs GMSH mesh generation via prepare_simulation(), reports element count
directly - WITHOUT launching Palace/MPI. Seconds-to-minutes, not the
~1hr+ this session's strap-at-2L/3 run took before dying with kappa=5071
and 48.7M elements.

Mirrors the eigenmode config in build_and_eigenmode_qw.py exactly (same
fine_mesh_components(['resonator1']) + fine_mesh_in_rectangle(strap) -
this design has no TL/ports at all, unlike the half-wave one, so there
is nothing to split).

Usage - check a single position:
    python mesh_test_qw_local.py --ground-pos 3.9527 \
        --palace-bin /home/ruiz/repo/spack/opt/spack/linux-zen3/palace-develop-blpbtxe43ne2or3yo7gou6ycmdwdxuzl/bin/palace

Usage - sweep several positions in one call to compare element counts
(this is the useful mode - run it across the whole strap range to spot
WHERE the pathology starts, before committing any single one to the
cluster):
    python mesh_test_qw_local.py --ground-pos-sweep \
        --palace-bin /home/ruiz/repo/spack/opt/spack/linux-zen3/palace-develop-blpbtxe43ne2or3yo7gou6ycmdwdxuzl/bin/palace
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def parse_ground_pos(s):
    if s.lower() == "none":
        return None
    return float(s)


def run_one(ground_pos_mm, args, design_module):
    build_design = design_module.build_design
    TOTAL_LENGTH_MM = design_module.TOTAL_LENGTH_MM

    t0 = time.time()
    design = build_design(ground_pos_mm=ground_pos_mm,
                           total_length_mm=args.total_length_mm or TOTAL_LENGTH_MM,
                           ground_overlap_um=args.strap_overlap_um)
    design.rebuild()

    from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
    import matplotlib
    matplotlib.use("Agg", force=True)

    eigen_options = {
        "mesh_refinement": 0,
        "dielectric_material": "silicon",
        "starting_freq": 5e9,
        "number_of_freqs": 2,
        "solns_to_save": 1,
        "solver_order": 2,
        "solver_tol": 1.0e-8,
        "solver_maxits": 200,
        "fillet_resolution": 12,
        "palace_dir": args.palace_bin,
        "num_cpus": 1,
    }

    outdir = os.path.expanduser("~/mesh_test_qw_scratch")
    os.makedirs(outdir, exist_ok=True)
    label = "none" if ground_pos_mm is None else f"{ground_pos_mm}mm"

    eigen_sim = PALACE_Eigenmode_Simulation(
        name=f"meshtest_qw_{label}",
        metal_design=design,
        sim_parent_directory=outdir,
        mode="simPC",
        meshing="GMSH",
        user_options=eigen_options,
        create_files=True,
    )
    eigen_sim.add_metallic(1)
    eigen_sim.add_ground_plane()

    eigen_sim.fine_mesh_components(
        ["resonator1"], min_size=args.bg_min_size_um * 1e-6,
        max_size=args.bg_max_size_um * 1e-6,
        taper_dist_min=10e-6, metals_only=False,
    )

    if ground_pos_mm is not None:
        poly_table = design.qgeometry.tables["poly"]
        strap_rows = poly_table[poly_table["name"].str.startswith("ground_strap_")]
        if len(strap_rows) == 0:
            print(f"    WARNING: no ground_strap_* geometry found for "
                  f"ground_pos={label}")
        else:
            strap_geom = strap_rows.iloc[0]["geometry"]
            minx, miny, maxx, maxy = strap_geom.bounds
            safe_margin_mm = args.strap_margin_um / 1000
            x1, y1 = (minx - safe_margin_mm) * 1e-3, (miny - safe_margin_mm) * 1e-3
            x2, y2 = (maxx + safe_margin_mm) * 1e-3, (maxy + safe_margin_mm) * 1e-3
            eigen_sim.fine_mesh_in_rectangle(
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
                min_size=args.strap_min_size_um * 1e-6, max_size=args.bg_max_size_um * 1e-6,
                taper_dist_min=args.strap_taper_dist_min_um * 1e-6,
                    taper_dist_max=args.strap_taper_dist_max_um * 1e-6,
            )

    t1 = time.time()
    eigen_sim.prepare_simulation()
    mesh_time = time.time() - t1

    import gmsh
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    elem_types, elem_tags, _ = gmsh.model.mesh.getElements()
    total_elements = sum(len(t) for t in elem_tags)

    gmsh.write(os.path.join(outdir, f"meshtest_qw_{label}", f"mesh_{label}.msh"))

    print(f"    ground_pos={label:>10}  elements={total_elements:>10,}  "
          f"nodes={len(node_tags):>10,}  mesh_time={mesh_time:6.1f}s  "
          f"total_time={time.time()-t0:6.1f}s")

    return total_elements


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ground-pos", type=parse_ground_pos, default=None,
                    help="mm from the SHORT end (sg1), or 'none' for baseline. "
                         "Ignored if --ground-pos-sweep is given.")
    p.add_argument("--ground-pos-sweep", action="store_true",
                    help="Check a spread of positions across the whole strap "
                         "range in one call, to spot where element count "
                         "blows up before committing a single position to "
                         "the cluster.")
    p.add_argument("--palace-bin", required=True,
                    help="Local Palace binary path. Never actually invoked "
                         "(no .run() call) but required at construction time.")
    p.add_argument("--total-length-mm", type=float, default=None,
                    help="Override TOTAL_LENGTH_MM. Default uses the design "
                         "module's own value.")
    p.add_argument("--bg-min-size-um", type=float, default=8.0)
    p.add_argument("--bg-max-size-um", type=float, default=100.0)
    p.add_argument("--strap-overlap-um", type=float, default=None,
                    help="How far each strap end extends PAST the outer gap "
                         "edge into the ground plane. None uses the design "
                         "module default. Sweep this if sliver elements "
                         "persist at the strap corners on filleted sections - "
                         "too small leaves a near-coincident edge, and the "
                         "straight-rectangle-vs-arc mismatch means corners on "
                         "a curve need more overlap than on a straight run.")
    p.add_argument("--strap-min-size-um", type=float, default=2.0,
                    help="Mesh size (um) INSIDE the small box around the "
                         "grounding strap only. 2um gives ~5 elements across "
                         "a 10um strap. Was previously hardcoded.")
    p.add_argument("--strap-taper-dist-min-um", type=float, default=2.0,
                    help="Distance (um) from the strap box within which the "
                         "full --strap-min-size-um applies.")
    p.add_argument("--strap-taper-dist-max-um", type=float, default=6.0,
                    help="Distance (um) beyond which the strap field relaxes "
                         "fully to --bg-max-size-um. Note the field returns "
                         "that max as a VALUE beyond this distance, so it "
                         "must taper to the background MAX or it caps the "
                         "whole domain - see the fine_mesh_in_rectangle call.")
    p.add_argument("--strap-margin-um", type=float, default=3.0)
    args = p.parse_args()

    import build_design_quarterwave as design_module

    if args.ground_pos_sweep:
        L = args.total_length_mm or design_module.TOTAL_LENGTH_MM
        margin = design_module.EDGE_MARGIN_MM
        # Deliberately includes 2L/3 explicitly - the exact position that
        # produced kappa=5071 / 48.7M elements on the cluster.
        positions = [None] + sorted(set([
            round(margin, 4), round(L / 6, 4), round(L / 3, 4),
            round(L / 2, 4), round(2 * L / 3, 4), round(5 * L / 6, 4),
            round(L - margin, 4),
        ]))
        print(f"Sweeping {len(positions)} positions for L={L}mm "
              f"(mesh {args.bg_min_size_um}/{args.bg_max_size_um}um)...\n")
        results = {}
        for pos in positions:
            try:
                results[pos] = run_one(pos, args, design_module)
            except Exception as e:
                print(f"    ground_pos={pos}  FAILED: {e}")
                results[pos] = None

        print("\n=== Summary ===")
        baseline = results.get(None)
        for pos, n in results.items():
            label = "none" if pos is None else f"{pos}mm"
            if n is None:
                print(f"    {label:>10}: FAILED")
            elif baseline:
                print(f"    {label:>10}: {n:>10,} elements  "
                      f"({n/baseline:5.1f}x baseline)")
            else:
                print(f"    {label:>10}: {n:>10,} elements")
    else:
        run_one(args.ground_pos, args, design_module)


if __name__ == "__main__":
    main()
