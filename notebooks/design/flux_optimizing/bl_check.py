#!/usr/bin/env python
"""Can a boundary-layer mesh (very fine ONLY at the thin PEC array edges, coarse
elsewhere) converge the aperture flux at a feasible tet count? The uniform-mesh
check showed the flux is genuinely field-under-resolved near the 2-8 um PEC
strips, not an integration artifact.

Each SGE task = one (point, strategy). Strategies:
  ref05  : uniform 0.5 um ribbon           -> reference "truth" (heavy)
  ref08  : uniform 0.8 um ribbon
  bl_a   : ribbon 2.0, array-edge 0.6, taper 2
  bl_b   : ribbon 2.0, array-edge 1.0, taper 3
  bl_c   : ribbon 1.5, array-edge 0.8, taper 1.5
  base3  : ribbon 3.0, array-edge 4.5  (current sweep)
"""
import argparse
import csv
import json
import os
import shutil
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_flux  # noqa: E402
from build_design import build_design, loop_aperture_polys  # noqa: E402

PHI0 = 2.067833848e-15

POINTS = {
    'nominal':       (35., 6., 25., 15., 4.),
    'thin_arm':      (35., 6., 25., 15., 2.),
    'small_all':     (30., 0., 14., 8., 2.),
    'big_loop_thin': (30., 0., 40., 24., 2.),
}
# (ribbon_min_um, array_edge_min_um, array_edge_taper_um)   None edge -> skip that call
STRATS = {
    'ref05': (0.5, None, None),
    'ref08': (0.8, None, None),
    'bl_a':  (2.0, 0.6, 2.0),
    'bl_b':  (2.0, 1.0, 3.0),
    'bl_c':  (1.5, 0.8, 1.5),
    'base3': (3.0, 4.5, 6.0),
}
TASKS = [(p, s) for p in POINTS for s in STRATS]


def one(name, strat, palace_bin, ncpu, workdir):
    ribbon, edge, taper = STRATS[strat]
    params = dict(zip(sweep_flux.AXIS_ORDER, POINTS[name]))
    sim_dir = os.path.join(workdir, 'flux')
    shutil.rmtree(sim_dir, ignore_errors=True)
    from SQDMetal.PALACE.Inductance_Simulation import PALACE_Inductance_Simulation
    t0 = time.time()

    d = build_design(chip_margin_um=250, chip_size_z_um=200, include_array_metal=True,
                     ground_near_end=True, **params)
    d.rebuild()
    ia = loop_aperture_polys(squid_gap_um=params['squid_gap_um'],
                             arm_semi_length_um=params['arm_semi_length_um'],
                             arm_width_um=params['arm_width_um'],
                             segment_length_um=params['segment_length_um'], grow_um=0.3)
    opts = {'mesh_refinement': 0, 'dielectric_material': 'silicon', 'solver_order': 1,
            'solver_tol': 1e-3, 'solver_maxits': 2000, 'fillet_resolution': 8,
            'solns_to_save': 0, 'palace_dir': palace_bin, 'num_cpus': ncpu}
    mag = PALACE_Inductance_Simulation(name='flux', metal_design=d,
                                       sim_parent_directory=workdir.rstrip('/') + '/',
                                       mode='simPC', meshing='GMSH', user_options=opts,
                                       create_files=True)
    mag.add_metallic(1)
    mag.add_ground_plane()
    mag.create_current_source_with_Uclip_on_Route('readout_1', 'start')
    mag.add_integration_area(ia)
    mag.set_xBoundary_as_proportion(0.05)
    mag.set_yBoundary_as_proportion(0.05)
    mag.set_zBoundary_as_proportion(z_prop_top=1.5, z_prop_bottom=0.0)
    x1, y1, x2, y2 = ia.bounds
    pad = 8e-6
    mag.fine_mesh_in_rectangle(x1 - pad, y1 - pad, x2 + pad, y2 + pad,
                               min_size=ribbon * 1e-6, max_size=60e-6, taper_dist_min=6e-6)
    if edge is not None:
        mag.fine_mesh_components(['array_1'], min_size=edge * 1e-6, max_size=60e-6,
                                 taper_dist_min=taper * 1e-6, metals_only=True)
    mag.prepare_simulation()
    cfg = json.load(open(mag._sim_config))
    cfg['Solver']['Linear']['KSPType'] = 'GMRES'
    cfg['Solver']['Linear']['MaxSize'] = 400
    json.dump(cfg, open(mag._sim_config, 'w'), indent=2)
    try:
        import gmsh
        ntet = len(gmsh.model.mesh.getElementsByType(4)[0])
    except Exception:
        ntet = None
    mag.run()
    raw = PALACE_Inductance_Simulation.retrieve_magnetostatic_data_from_file(mag._output_data_dir)
    phi = float(np.ravel(raw['surface_Flux'])[0])
    i_inc = float(np.ravel(raw['terminal_I'])[0])
    shutil.rmtree(sim_dir, ignore_errors=True)
    return {'point': name, 'strat': strat, 'Phi0_per_mA': (phi / i_inc) / PHI0 * 1e-3,
            'n_tet': ntet, 'seconds': round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--index', type=int, required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--palace-bin', required=True)
    ap.add_argument('--num-cpus', type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    sweep_flux._patch_gmsh_multipoly()
    name, strat = TASKS[a.index - 1]
    wd = os.path.join(a.outdir, f'work_{a.index:02d}')
    os.makedirs(wd, exist_ok=True)
    print(f"[{a.index}] {name} / {strat}", flush=True)
    try:
        row = one(name, strat, a.palace_bin, a.num_cpus, wd)
    except Exception as e:
        row = {'point': name, 'strat': strat, 'Phi0_per_mA': '', 'n_tet': '',
               'seconds': '', 'error': f"{type(e).__name__}: {e}"[:200]}
    with open(os.path.join(a.outdir, f'row_{a.index:02d}.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    print(f"[{a.index}] {row}", flush=True)
    shutil.rmtree(wd, ignore_errors=True)


if __name__ == '__main__':
    print(f"{len(TASKS)} tasks", flush=True)
    main()
