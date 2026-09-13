#!/usr/bin/env python
"""Magnetostatic self-inductance-per-unit-length of a straight CPW segment: a
current is driven at the open end (via a U-clip that locally shorts it to
ground, for the sim only) and returns through ground via the real short at the
far end, so the closed loop it forms encloses the CPW gap along the trace's
whole length. Flux through that loop, divided by the driven current, is the
trace's self-inductance; dividing by the (margin-trimmed) length gives L'.

    python run_inductance.py --outdir out_ind --palace-bin $PALACE_BIN \
        --cpw-width-um 20 --cpw-gap-um 12 --num-cpus 8

Cluster job only -- run via submit_inductance.sh. Do not run this locally
(PALACE/FEM solves always go to the CNEA cluster).

Mesh sizing: base fine-mesh element size is HALF the smallest lateral feature
(min(cpw_width, cpw_gap)) -- see run_capacitance.py's docstring for the same
rule and why, for a FULL3D thickness sweep, the through-thickness floor
follows scripts/python/island_cap/run_cap_island.py's convention instead of
literally halving a 20-300nm thickness. AMR (enable_mesh_refinement) then
refines from there toward a target relative error.
"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import (build_design, integration_area_gap_strips,
                          OTG_NAME, STG_NAME, ROUTE_NAME, LENGTH_UM)


def _patch_gmsh_multipoly():
    """SQDMetal's GMSH builder chokes on a MultiPolygon input ("'MultiPolygon'
    object has no attribute 'exterior'") -- hit here on the metal union (the
    trace + terminations + U-clip don't always merge into a single Polygon).
    Same fix as notebooks/design/flux_optimizing/run_flux.py: flatten any
    MultiPolygon into its constituent Polygons before GMSH sees it.
    """
    from SQDMetal.PALACE.Utilities.GMSH_Geometry_Builder import GMSH_Geometry_Builder
    orig = GMSH_Geometry_Builder._create_gmsh_geometry_from_shapely_polygons

    def patched(self, polygons):
        flat = []
        for p in polygons:
            flat.extend(list(p.geoms) if p.geom_type == 'MultiPolygon' else [p])
        return orig(self, flat)

    GMSH_Geometry_Builder._create_gmsh_geometry_from_shapely_polygons = patched


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--palace-bin', required=True)
    ap.add_argument('--cpw-width-um', type=float, required=True)
    ap.add_argument('--cpw-gap-um', type=float, required=True)
    ap.add_argument('--length-um', type=float, default=LENGTH_UM)
    ap.add_argument('--margin-um', type=float, default=None,
                    help='trimmed off each x-end of the flux integral; default '
                         'max(5*cpw_gap, 30um), see build_design.integration_area_gap_strips')
    ap.add_argument('--num-cpus', type=int, default=8)
    ap.add_argument('--loop-mesh-um', type=float, default=None,
                    help='fine-mesh target across the gap; default 0.5*min(width,gap)')
    ap.add_argument('--full3d', action='store_true',
                    help='real metal thickness instead of the default zero-thickness sheet')
    ap.add_argument('--metal-thickness-nm', type=float, default=100.0,
                    help='only used with --full3d')
    ap.add_argument('--solver-order', type=int, default=1)
    ap.add_argument('--solver-tol', type=float, default=1e-8)
    ap.add_argument('--solver-maxits', type=int, default=500)
    ap.add_argument('--ksp', default='CG')
    ap.add_argument('--amr', dest='amr', action='store_true', default=True)
    ap.add_argument('--no-amr', dest='amr', action='store_false')
    ap.add_argument('--amr-iters', type=int, default=4)
    ap.add_argument('--amr-tol', type=float, default=1e-2)
    ap.add_argument('--amr-max-dofs', type=int, default=3_000_000)
    ap.add_argument('--amr-dorfler', type=float, default=0.6)
    ap.add_argument('--tag', default=None)
    args = ap.parse_args()

    tag = args.tag or f"w{args.cpw_width_um:g}_g{args.cpw_gap_um:g}"
    os.makedirs(args.outdir, exist_ok=True)
    _patch_gmsh_multipoly()

    t0 = time.time()
    from qiskit_metal import designs  # noqa: F401
    import matplotlib
    matplotlib.use("Agg", force=True)
    from SQDMetal.PALACE.Inductance_Simulation import PALACE_Inductance_Simulation

    # far_end='short': a REAL galvanic short, needed to close the current loop
    # (see build_design.build_design's docstring -- do NOT change this to
    # 'open', that's the capacitance sim's geometry, not this one's).
    design = build_design(args.cpw_width_um, args.cpw_gap_um, length_um=args.length_um,
                          far_end='short')
    print(f"[{tag}] geometry built  chip={design.chips.main.size.size_x} x "
          f"{design.chips.main.size.size_y}  FULL3D={args.full3d} "
          f"thickness={args.metal_thickness_nm if args.full3d else 'n/a'}nm", flush=True)

    int_area = integration_area_gap_strips(args.cpw_width_um, args.cpw_gap_um,
                                           length_um=args.length_um,
                                           margin_um=args.margin_um)
    margin_um = args.margin_um if args.margin_um is not None else max(5.0 * args.cpw_gap_um, 30.0)
    integration_length_um = args.length_um - 2 * margin_um

    opts = {
        'mesh_refinement': 0,
        'dielectric_material': 'silicon',
        'solver_order': args.solver_order,
        'solver_tol': args.solver_tol,
        'solver_maxits': args.solver_maxits,
        'fillet_resolution': 8,
        'palace_dir': args.palace_bin,
        'num_cpus': args.num_cpus,
    }
    mag = PALACE_Inductance_Simulation(
        name=f'ci_{tag}', metal_design=design, sim_parent_directory=args.outdir.rstrip('/') + '/',
        mode='simPC', meshing='GMSH', user_options=opts, create_files=True)

    mag.add_metallic(1)
    mag.add_ground_plane()
    if args.full3d:
        mag.enforce_full_3D_simulation(metal_thickness=args.metal_thickness_nm * 1e-9,
                                       substrate_trenching=0.0)
    mag.create_current_source_with_Uclip_on_Route(ROUTE_NAME, 'start')
    mag.add_integration_area(int_area)

    mag.set_xBoundary_as_proportion(0.05)
    mag.set_yBoundary_as_proportion(0.05)
    mag.set_zBoundary_as_proportion(z_prop_top=1.5, z_prop_bottom=0.0)

    # half the smallest lateral feature, floor'd against the metal thickness
    # for FULL3D (see module docstring) -- terminations MUST be in the
    # fine-mesh list (sqdmetal-termination-gap-meshing lesson) or the open
    # end's isolating notch collapses and the loop this sim depends on never
    # forms.
    lateral_min_um = 0.5 * min(args.cpw_width_um, args.cpw_gap_um)
    loop_mesh_um = args.loop_mesh_um or max(
        lateral_min_um, args.metal_thickness_nm * 1e-3 if args.full3d else 0.0, 0.15)
    x1, y1, x2, y2 = int_area.bounds
    pad = 3.0 * args.cpw_gap_um * 1e-6
    mag.fine_mesh_in_rectangle(x1 - pad, y1 - pad, x2 + pad, y2 + pad,
                               min_size=loop_mesh_um * 1e-6, max_size=80e-6,
                               taper_dist_min=6 * args.cpw_gap_um * 1e-6)
    mag.fine_mesh_components([ROUTE_NAME, OTG_NAME, STG_NAME],
                             min_size=loop_mesh_um * 1e-6, max_size=80e-6,
                             taper_dist_min=6e-6, metals_only=False)
    print(f"mesh: lateral_min={lateral_min_um:.3f}um  loop_mesh={loop_mesh_um:.3f}um", flush=True)

    if args.amr:
        mag.enable_mesh_refinement(
            num_iterations=args.amr_iters, max_DoFs=args.amr_max_dofs, tolerance=args.amr_tol,
            Dorfler_marking_fraction=args.amr_dorfler,
            save_iterations_data=True, save_iterations_mesh=False, nonconformal=False)
        print(f"AMR enabled: iters={args.amr_iters} tol={args.amr_tol} "
              f"max_DoFs={args.amr_max_dofs:,} Dorfler={args.amr_dorfler}", flush=True)

    mag.prepare_simulation()

    # match run_flux.py: SQDMetal hardcodes CG, which can stall on topologically
    # nontrivial geometry. A single straight loop has no such handle, so CG should
    # be fine, but keep the override available (--ksp GMRES) if a run stalls.
    if args.ksp != 'CG':
        cfg_path = mag._sim_config
        cfg = json.load(open(cfg_path))
        cfg['Solver']['Linear']['KSPType'] = args.ksp
        json.dump(cfg, open(cfg_path, 'w'), indent=2)
    print(f"[solver] KSP={args.ksp}", flush=True)

    try:
        import gmsh
        n_tet = len(gmsh.model.mesh.getElementsByType(4)[0])
        n_nod = len(gmsh.model.mesh.getNodes()[0])
        print(f"[mesh] tets={n_tet:,}  nodes={n_nod:,}", flush=True)
    except Exception as e:
        n_tet = n_nod = None
        print(f"[mesh] count failed: {e}", flush=True)

    t_solve = time.time()
    mag.run()
    solve_s = time.time() - t_solve

    def _L_total(output_dir):
        raw = PALACE_Inductance_Simulation.retrieve_magnetostatic_data_from_file(output_dir)
        phi = float(np.ravel(raw['surface_Flux'])[0])
        i_inc = float(np.ravel(raw['terminal_I'])[0])
        return phi / i_inc, phi, i_inc

    # per-AMR-iteration convergence of L itself (not just Palace's internal FEM
    # error estimator) -- same pattern as squid_array_coupling/run_amr.py.
    convergence = []
    if args.amr:
        for it_dir in sorted(glob.glob(os.path.join(mag._output_data_dir, "iteration*"))):
            try:
                L_it, _, _ = _L_total(it_dir)
                convergence.append({"iteration": os.path.basename(it_dir), "L_nH": L_it * 1e9})
            except Exception as e:
                convergence.append({"iteration": os.path.basename(it_dir), "error": str(e)})
        for i in range(1, len(convergence)):
            if "L_nH" in convergence[i] and "L_nH" in convergence[i - 1]:
                convergence[i]["rel_change_vs_prev"] = (
                    abs(convergence[i]["L_nH"] - convergence[i - 1]["L_nH"]) / convergence[i]["L_nH"])
        print("AMR convergence (L_nH per iteration): " + json.dumps(convergence, indent=2), flush=True)

    L_total, phi, i_inc = _L_total(mag._output_data_dir)
    Lprime_H_per_m = L_total / (integration_length_um * 1e-6)
    Lprime_nH_per_um = Lprime_H_per_m * 1e3                  # 1 H/m = 1e3 nH/um

    result = {
        'tag': tag, 'cpw_width_um': args.cpw_width_um, 'cpw_gap_um': args.cpw_gap_um,
        'length_um': args.length_um, 'margin_um': margin_um,
        'integration_length_um': integration_length_um,
        'full3d': args.full3d,
        'metal_thickness_nm': args.metal_thickness_nm if args.full3d else None,
        'Phi_Wb': phi, 'I_A': i_inc, 'L_total_H': L_total,
        'Lprime_H_per_m': Lprime_H_per_m, 'Lprime_nH_per_um': Lprime_nH_per_um,
        'n_tet': n_tet, 'n_nod': n_nod,
        'amr': {'enabled': args.amr, 'iterations_requested': args.amr_iters,
               'tolerance': args.amr_tol, 'max_dofs': args.amr_max_dofs,
               'convergence': convergence} if args.amr else {'enabled': False},
        'solve_seconds': round(solve_s, 1), 'total_seconds': round(time.time() - t0, 1),
    }
    os.makedirs(os.path.join(args.outdir, 'results'), exist_ok=True)
    with open(os.path.join(args.outdir, 'results', f'{tag}.json'), 'w') as f:
        json.dump(result, f, indent=2)
    print("[result] " + json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
