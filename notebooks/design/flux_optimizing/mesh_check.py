#!/usr/bin/env python
"""Mesh-convergence check for the flux sweep: does a flat sweep-grade mesh
(LOOP_MESH_UM = 3) bias the flux for small-feature geometries, changing the
ranking? Runs a handful of grid points, each at several mesh levels, one
(point, mesh) per SGE array task.

    python mesh_check.py --index $SGE_TASK_ID --outdir meshchk \
        --palace-bin ... --num-cpus $NSLOTS
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sweep_flux  # noqa: E402

# (route_clearance, segment_length, squid_gap, arm_semi_length, arm_width)  [um]
POINTS = {
    'nominal':        (35., 6., 25., 15., 4.),
    'thin_arm':       (35., 6., 25., 15., 2.),
    'thick_arm':      (35., 6., 25., 15., 8.),
    'small_all':      (30., 0., 14., 8., 2.),
    'small_loop_thin':(35., 6., 14., 8., 2.),
    'big_loop_thin':  (30., 0., 40., 24., 2.),
    'big_all':        (40., 12., 40., 24., 8.),
    'mid':            (30., 12., 25., 24., 4.),
}
MESHES = [3.0, 2.0, 1.5]
FINE_EXTRA = {'thin_arm', 'small_all', 'small_loop_thin', 'big_loop_thin'}  # + 1.0 um

# build the flat task list: (name, mesh)
TASKS = []
for name in POINTS:
    for m in MESHES:
        TASKS.append((name, m))
    if name in FINE_EXTRA:
        TASKS.append((name, 1.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--index', type=int, required=True)   # 1..len(TASKS)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--palace-bin', required=True)
    ap.add_argument('--num-cpus', type=int, default=8)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    sweep_flux._patch_gmsh_multipoly()

    name, mesh = TASKS[args.index - 1]
    params = dict(zip(sweep_flux.AXIS_ORDER, POINTS[name]))
    sweep_flux.LOOP_MESH_UM = mesh            # override the module constant
    workdir = os.path.join(args.outdir, f'work_{args.index:02d}')
    os.makedirs(workdir, exist_ok=True)

    print(f"[{args.index}] point={name} {params}  mesh={mesh} um", flush=True)
    row = sweep_flux.run_one(params, args.palace_bin, args.num_cpus, workdir)
    row = {'point': name, 'mesh_um': mesh, **row}

    csv_path = os.path.join(args.outdir, f'row_{args.index:02d}.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(row))
        w.writeheader()
        w.writerow(row)
    print(f"[{args.index}] -> Phi0/mA={row.get('Phi0_per_mA')} n_tet={row.get('n_tet')} "
          f"n_ksp={row.get('n_ksp')} err={row.get('error')}", flush=True)

    import shutil
    shutil.rmtree(workdir, ignore_errors=True)


if __name__ == '__main__':
    print(f"{len(TASKS)} tasks total", flush=True)
    main()
