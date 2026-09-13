#!/usr/bin/env python
"""Magnetostatic flux-per-amp of the readout hairpin through the SQUID-loop
apertures, with the array metal either IN or OUT of the EM model.

Mirrors inductance_sim.ipynb. The point: does dropping the array metal (which
makes the loops topological handles that stall the AMS solver) change the flux?

    python run_flux.py --outdir out_array            # array metal in
    python run_flux.py --outdir out_noarray --no-array
"""
import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import build_design, loop_aperture_polys  # noqa: E402

PHI0 = 2.067833848e-15  # Wb


def _patch_gmsh_multipoly():
    """SQDMetal's GMSH builder chokes on a MultiPolygon dielectric-gap piece
    (which our 23-aperture integration area produces). Flatten MultiPolygons
    into their constituent Polygons -- exactly what the surrounding code already
    does for 'left-over' pieces."""
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
    ap.add_argument('--num-cpus', type=int, default=32)
    ap.add_argument('--no-array', action='store_true',
                    help='leave the SQUID-array metal out of the EM model')
    ap.add_argument('--loop-mesh-um', type=float, default=2.5)
    ap.add_argument('--ribbon-taper-um', type=float, default=6.0)
    ap.add_argument('--ribbon-pad-um', type=float, default=8.0)
    ap.add_argument('--solver-order', type=int, default=1)
    ap.add_argument('--solver-tol', type=float, default=1e-6)
    ap.add_argument('--solver-maxits', type=int, default=2000)
    ap.add_argument('--ksp', default='GMRES')
    ap.add_argument('--chip-margin-um', type=float, default=250.0)
    ap.add_argument('--chip-z-um', type=float, default=200.0)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    _patch_gmsh_multipoly()

    from SQDMetal.PALACE.Inductance_Simulation import PALACE_Inductance_Simulation

    with_array = not args.no_array
    t0 = time.time()

    design = build_design(chip_margin_um=args.chip_margin_um,
                          chip_size_z_um=args.chip_z_um,
                          include_array_metal=with_array)
    design.rebuild()
    sz = design.chips.main.size
    print(f"[build] array_metal={with_array}  components={list(design.components)}  "
          f"chip={sz.size_x} x {sz.size_y} x {sz.size_z}", flush=True)

    # identical integration area either way: the 23 loop apertures, from formula
    int_area = loop_aperture_polys()

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
        name='flux', metal_design=design, sim_parent_directory=args.outdir.rstrip('/') + '/',
        mode='simPC', meshing='GMSH', user_options=opts, create_files=True)

    mag.add_metallic(1)
    mag.add_ground_plane()
    mag.create_current_source_with_Uclip_on_Route('readout_1', 'start')
    mag.add_integration_area(int_area)

    mag.set_xBoundary_as_proportion(0.05)
    mag.set_yBoundary_as_proportion(0.05)
    mag.set_zBoundary_as_proportion(z_prop_top=1.5, z_prop_bottom=0.0)

    x1, y1, x2, y2 = int_area.bounds
    pad = args.ribbon_pad_um * 1e-6
    mag.fine_mesh_in_rectangle(x1 - pad, y1 - pad, x2 + pad, y2 + pad,
                               min_size=args.loop_mesh_um * 1e-6, max_size=60e-6,
                               taper_dist_min=args.ribbon_taper_um * 1e-6)
    if with_array:
        mag.fine_mesh_components(['array_1'], min_size=1.5 * args.loop_mesh_um * 1e-6,
                                 max_size=60e-6, taper_dist_min=6e-6, metals_only=True)

    mag.prepare_simulation()

    # SQDMetal hardcodes CG; it stalls with the array present -> GMRES
    cfg_path = mag._sim_config
    cfg = json.load(open(cfg_path))
    cfg['Solver']['Linear']['KSPType'] = args.ksp
    json.dump(cfg, open(cfg_path, 'w'), indent=2)
    print(f"[solver] {cfg['Solver']['Linear']}", flush=True)

    try:
        import gmsh
        n_tet = len(gmsh.model.mesh.getElementsByType(4)[0])
        n_nod = len(gmsh.model.mesh.getNodes()[0])
        print(f"[mesh] tets={n_tet:,}  nodes={n_nod:,}", flush=True)
    except Exception as e:
        n_tet = n_nod = None
        print(f"[mesh] count failed: {e}", flush=True)

    t_solve = time.time()
    flux_per_amp = np.array(mag.run()).ravel()
    solve_s = time.time() - t_solve

    log = os.path.join(mag._output_data_dir, 'out.log')
    n_iters = None
    converged = None
    if os.path.exists(log):
        txt = open(log).read()
        ks = [l for l in txt.splitlines() if 'KSP residual' in l]
        n_iters = len(ks)
        converged = ('solver converged' in txt.lower()) and ('did not converge' not in txt.lower())

    raw = PALACE_Inductance_Simulation.retrieve_magnetostatic_data_from_file(mag._output_data_dir)
    phi = float(np.ravel(raw['surface_Flux'])[0])
    i_inc = float(np.ravel(raw['terminal_I'])[0])

    result = {
        'with_array_metal': with_array,
        'Phi_mag_Wb': phi,
        'I_inc_A': i_inc,
        'flux_per_amp_Wb_per_A': phi / i_inc,
        'Phi0_per_mA': (phi / i_inc) / PHI0 * 1e-3,
        'n_tet': n_tet, 'n_nod': n_nod,
        'n_ksp_iters': n_iters, 'converged': converged,
        'solve_seconds': round(solve_s, 1),
        'total_seconds': round(time.time() - t0, 1),
        'solver': cfg['Solver']['Linear'],
    }
    with open(os.path.join(args.outdir, 'result.json'), 'w') as f:
        json.dump(result, f, indent=2)
    print("[result] " + json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
