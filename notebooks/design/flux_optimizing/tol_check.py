#!/usr/bin/env python
"""Is the flux non-convergence a SOLVER issue (the 23 open-loop cohomology modes
converge at unpreconditioned GMRES rate, ~0.985/iter) rather than a mesh issue?
Test: same geometry at 2 meshes x 2 tolerances, huge maxits. If the tight-tol
values agree across mesh, tighter tol + more iters is the (expensive) fix.
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
POINTS = {'nominal': (35., 6., 25., 15., 4.), 'small_all': (30., 0., 14., 8., 2.)}
MESHES = [3.0, 1.5]
TOLS = [1e-3, 1e-7]
TASKS = [(p, m, t) for p in POINTS for m in MESHES for t in TOLS]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--index', type=int, required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--palace-bin', required=True)
    ap.add_argument('--num-cpus', type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    sweep_flux._patch_gmsh_multipoly()
    from SQDMetal.PALACE.Inductance_Simulation import PALACE_Inductance_Simulation

    name, mesh, tol = TASKS[a.index - 1]
    params = dict(zip(sweep_flux.AXIS_ORDER, POINTS[name]))
    wd = os.path.join(a.outdir, f'work_{a.index:02d}')
    os.makedirs(wd, exist_ok=True)
    sim_dir = os.path.join(wd, 'flux')
    t0 = time.time()
    row = {'point': name, 'mesh_um': mesh, 'tol': tol}
    try:
        d = build_design(chip_margin_um=250, chip_size_z_um=200, include_array_metal=True,
                         ground_near_end=True, **params)
        d.rebuild()
        ia = loop_aperture_polys(squid_gap_um=params['squid_gap_um'],
                                 arm_semi_length_um=params['arm_semi_length_um'],
                                 arm_width_um=params['arm_width_um'],
                                 segment_length_um=params['segment_length_um'], grow_um=0.3)
        opts = {'mesh_refinement': 0, 'dielectric_material': 'silicon', 'solver_order': 1,
                'solver_tol': tol, 'solver_maxits': 12000, 'fillet_resolution': 8,
                'solns_to_save': 0, 'palace_dir': a.palace_bin, 'num_cpus': a.num_cpus}
        mag = PALACE_Inductance_Simulation(name='flux', metal_design=d,
                                           sim_parent_directory=wd.rstrip('/') + '/',
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
                                   min_size=mesh * 1e-6, max_size=60e-6, taper_dist_min=6e-6)
        mag.fine_mesh_components(['array_1'], min_size=1.5 * mesh * 1e-6, max_size=60e-6,
                                 taper_dist_min=6e-6, metals_only=True)
        mag.prepare_simulation()
        cfg = json.load(open(mag._sim_config))
        cfg['Solver']['Linear']['KSPType'] = 'GMRES'
        cfg['Solver']['Linear']['MaxSize'] = 1000
        json.dump(cfg, open(mag._sim_config, 'w'), indent=2)
        try:
            import gmsh
            row['n_tet'] = len(gmsh.model.mesh.getElementsByType(4)[0])
        except Exception:
            row['n_tet'] = ''
        mag.run()
        raw = PALACE_Inductance_Simulation.retrieve_magnetostatic_data_from_file(mag._output_data_dir)
        phi = float(np.ravel(raw['surface_Flux'])[0])
        i_inc = float(np.ravel(raw['terminal_I'])[0])
        row['Phi0_per_mA'] = phi / i_inc / PHI0 * 1e-3
        txt = open(os.path.join(mag._output_data_dir, 'out.log')).read()
        row['n_ksp'] = sum('KSP residual' in ln for ln in txt.splitlines())
        row['converged'] = 'solver converged' in txt.lower()
        # last residual
        res = [ln for ln in txt.splitlines() if 'KSP residual' in ln]
        row['last_resid'] = res[-1].split('=')[-1].strip() if res else ''
    except Exception as e:
        row['error'] = f"{type(e).__name__}: {e}"[:200]
    row['seconds'] = round(time.time() - t0, 1)
    shutil.rmtree(sim_dir, ignore_errors=True)
    with open(os.path.join(a.outdir, f'row_{a.index:02d}.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    print(row, flush=True)
    shutil.rmtree(wd, ignore_errors=True)


if __name__ == '__main__':
    print(f"{len(TASKS)} tasks", flush=True)
    main()
