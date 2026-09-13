#!/usr/bin/env python
"""Does moving the flux-integration surface off the singular metal edge fix the
mesh non-convergence? Sweeps (point, mesh, aperture-grow); one per SGE task.

grow modes:
  edge      : +0.3 um  (current sweep -- boundary just past the inner metal edge,
              in the B-field edge singularity)
  mid       : +1.5 um
  center    : +arm_width/2  (boundary at the wire centreline, where B_z ~ 0)
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_flux  # noqa: E402

POINTS = {
    'nominal':       (35., 6., 25., 15., 4.),
    'thin_arm':      (35., 6., 25., 15., 2.),
    'thick_arm':     (35., 6., 25., 15., 8.),
    'small_all':     (30., 0., 14., 8., 2.),
    'big_loop_thin': (30., 0., 40., 24., 2.),
    'mid':           (30., 12., 25., 24., 4.),
}
MESHES = [3.0, 2.0]
GROWS = ['edge', 'mid', 'center']

TASKS = [(p, m, g) for p in POINTS for m in MESHES for g in GROWS]


def grow_um_for(mode, params):
    if mode == 'edge':
        return 0.3
    if mode == 'mid':
        return 1.5
    if mode == 'center':
        return 0.5 * params['arm_width_um']
    raise ValueError(mode)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--index', type=int, required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--palace-bin', required=True)
    ap.add_argument('--num-cpus', type=int, default=8)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sweep_flux._patch_gmsh_multipoly()

    name, mesh, growmode = TASKS[args.index - 1]
    params = dict(zip(sweep_flux.AXIS_ORDER, POINTS[name]))
    grow = grow_um_for(growmode, params)

    sweep_flux.LOOP_MESH_UM = mesh
    sweep_flux.GROW_TRIES = (grow, grow + 1.0, grow + 2.5)   # fixed base grow, retry adds

    workdir = os.path.join(args.outdir, f'work_{args.index:02d}')
    os.makedirs(workdir, exist_ok=True)
    print(f"[{args.index}] {name} {params} mesh={mesh} grow={growmode}({grow:.2f}um)", flush=True)

    row = sweep_flux.run_one(params, args.palace_bin, args.num_cpus, workdir)
    row = {'point': name, 'mesh_um': mesh, 'grow_mode': growmode, **row}

    with open(os.path.join(args.outdir, f'row_{args.index:02d}.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    print(f"[{args.index}] -> Phi0/mA={row.get('Phi0_per_mA')} n_tet={row.get('n_tet')} "
          f"err={row.get('error')}", flush=True)

    import shutil
    shutil.rmtree(workdir, ignore_errors=True)


if __name__ == '__main__':
    print(f"{len(TASKS)} tasks", flush=True)
    main()
