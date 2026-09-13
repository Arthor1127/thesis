#!/usr/bin/env python
"""One shard of the 5-D flux-threading sweep.

Grid axes (all um):  route_clearance, segment_length, squid_gap,
                     arm_semi_length, arm_width  ->  4*5*6*5*4 = 2400 points.

Each point: build the array + hairpin design (array metal IN, array near-end tied
to ground), magnetostatic Palace solve (GMRES), record Sum_i Phi_i / I.
Mesh/paraview deleted after each point. Rows appended (and flushed) to
<outdir>/shard_<K>.csv as they finish, so a killed shard resumes where it
stopped on resubmit.

    python sweep_flux.py --shard 1 --nshards 200 --outdir sweep \
        --palace-bin /path/to/palace --num-cpus 8
"""
import argparse
import csv
import gc
import itertools
import os
import shutil
import sys
import time
import traceback

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

PHI0 = 2.067833848e-15  # Wb

# --- the grid -------------------------------------------------------------
# Trimmed to ~2400 pts to fit a ~9-12 h budget at the converged boundary-layer
# mesh (see below). arm_width starts at 3 um (2 um doesn't mesh-converge and is
# the lowest-flux corner anyway).
AXES = {
    'route_clearance_um': np.array([30., 33., 36., 40.]),
    'segment_length_um':  np.array([0., 3., 6., 9., 12.]),
    'squid_gap_um':       np.round(np.linspace(14., 40., 6), 1),
    'arm_semi_length_um': np.array([8., 12., 16., 20., 24.]),
    'arm_width_um':       np.array([3., 4.5, 6., 8.]),
}   # 4 * 5 * 6 * 5 * 4 = 2400
AXIS_ORDER = list(AXES)
FIELDS = AXIS_ORDER + ['flux_per_amp', 'Phi0_per_mA', 'Phi_mag_Wb', 'I_inc_A',
                       'grow_um', 'n_tet', 'n_ksp', 'converged', 'seconds', 'error']

# integration-aperture outward-grow tries: exact metal-edge coincidence makes GMSH
# fail ("Could not create line"); a larger grow moves the surface clear of it.
GROW_TRIES = (0.3, 1.0, 2.5)

# --- solver / mesh ------------------------------------------------------------
# Uniform 1.5 um ribbon over the loop band + array metal. From the convergence
# study: nominal 1.5 um = 13.71 Phi0/mA vs converged ~13.95 (-1.7%, a near-flat
# offset -> ranking is clean), ~1.75 M tets, ~250 s. The "bl_b" boundary-layer
# mesh was converged to <1% but 30-60 min/point (its 1 um PEC-edge shell sharpens
# the open-loop cohomology modes AND adds DoF -> GMRES needs 2000-3000 iters).
# solver_tol 1e-3 == 1e-7 to 0.06% here. Points that hit SOLVER_MAXITS are still
# recorded (converged=False) -- flag and inspect, don't let them grind forever.
RIBBON_MESH_UM = 1.5
ARRAY_EDGE_UM = 1.5
ARRAY_EDGE_TAPER_UM = 6.0
SOLVER_TOL = 1e-3
SOLVER_MAXITS = 2000


def all_points():
    return list(itertools.product(*[AXES[a] for a in AXIS_ORDER]))


def _patch_gmsh_multipoly():
    from SQDMetal.PALACE.Utilities.GMSH_Geometry_Builder import GMSH_Geometry_Builder
    orig = GMSH_Geometry_Builder._create_gmsh_geometry_from_shapely_polygons

    def patched(self, polygons):
        flat = []
        for p in polygons:
            flat.extend(list(p.geoms) if p.geom_type == 'MultiPolygon' else [p])
        return orig(self, flat)
    GMSH_Geometry_Builder._create_gmsh_geometry_from_shapely_polygons = patched


def run_one(params, palace_bin, num_cpus, workdir, solver_maxits=SOLVER_MAXITS,
            gmres_maxsize=250):
    """One grid point -> dict of FIELDS. Never raises; errors go in 'error'.
    Retries with a larger integration-aperture grow if GMSH hits the
    coincident-edge failure."""
    row = {k: params[k] for k in AXIS_ORDER}
    row.update({f: '' for f in FIELDS if f not in row})
    t0 = time.time()
    retryable = ('Could not create line', 'surface-F.csv', 'No such file')
    for grow in GROW_TRIES:
        r = _attempt(params, palace_bin, num_cpus, workdir, grow,
                    solver_maxits=solver_maxits, gmres_maxsize=gmres_maxsize)
        r['seconds'] = round(time.time() - t0, 1)
        if r.get('flux_per_amp', '') != '':
            return r
        if not any(s in str(r.get('error', '')) for s in retryable):
            return r          # a different failure -- retrying grow won't help
        row = r               # keep the last error if all grows fail
    return row


def _attempt(params, palace_bin, num_cpus, workdir, grow_um,
            solver_maxits=SOLVER_MAXITS, gmres_maxsize=250):
    from build_design import build_design, loop_aperture_polys
    from SQDMetal.PALACE.Inductance_Simulation import PALACE_Inductance_Simulation

    row = {k: params[k] for k in AXIS_ORDER}
    row.update({f: '' for f in FIELDS if f not in row})
    row['grow_um'] = grow_um
    sim_dir = os.path.join(workdir, 'flux')
    try:
        if os.path.isdir(sim_dir):
            shutil.rmtree(sim_dir, ignore_errors=True)

        design = build_design(
            chip_margin_um=250, chip_size_z_um=200, include_array_metal=True,
            ground_near_end=True,   # array's near (elbow-side) end tied to ground
            array_ground_margin_um=None,  # THIS sweep keeps SquidArray's conformal
            #   buffer(ground_gap) pocket (with the between-loop ground tabs) for
            #   consistency across all 2400 points -- the switch to a clean
            #   rectangular pocket is evaluated separately in the corroboration runs.
            route_clearance_um=params['route_clearance_um'],
            segment_length_um=params['segment_length_um'],
            squid_gap_um=params['squid_gap_um'],
            arm_semi_length_um=params['arm_semi_length_um'],
            arm_width_um=params['arm_width_um'])
        design.rebuild()

        int_area = loop_aperture_polys(
            squid_gap_um=params['squid_gap_um'],
            arm_semi_length_um=params['arm_semi_length_um'],
            arm_width_um=params['arm_width_um'],
            segment_length_um=params['segment_length_um'],
            grow_um=grow_um)

        opts = {'mesh_refinement': 0, 'dielectric_material': 'silicon',
                'solver_order': 1, 'solver_tol': SOLVER_TOL,
                'solver_maxits': solver_maxits, 'fillet_resolution': 8,
                'solns_to_save': 0,   # no paraview VTU -- we only need the CSVs
                'palace_dir': palace_bin, 'num_cpus': num_cpus}
        mag = PALACE_Inductance_Simulation(
            name='flux', metal_design=design, sim_parent_directory=workdir.rstrip('/') + '/',
            mode='simPC', meshing='GMSH', user_options=opts, create_files=True)
        mag.add_metallic(1)
        mag.add_ground_plane()
        mag.create_current_source_with_Uclip_on_Route('readout_1', 'start')
        mag.add_integration_area(int_area)
        mag.set_xBoundary_as_proportion(0.05)
        mag.set_yBoundary_as_proportion(0.05)
        mag.set_zBoundary_as_proportion(z_prop_top=1.5, z_prop_bottom=0.0)
        x1, y1, x2, y2 = int_area.bounds
        pad = 8e-6
        mag.fine_mesh_in_rectangle(x1 - pad, y1 - pad, x2 + pad, y2 + pad,
                                   min_size=RIBBON_MESH_UM * 1e-6, max_size=60e-6,
                                   taper_dist_min=6e-6)
        mag.fine_mesh_components(['array_1'], min_size=ARRAY_EDGE_UM * 1e-6,
                                 max_size=60e-6,
                                 taper_dist_min=ARRAY_EDGE_TAPER_UM * 1e-6,
                                 metals_only=True)
        mag.prepare_simulation()

        import json
        cfg = json.load(open(mag._sim_config))
        cfg['Solver']['Linear']['KSPType'] = 'GMRES'
        # GMRES restart size. Bigger = fewer restarts on the floating-island
        # near-null-space (slow points need 600-1600 iters), but each vector is
        # n_DoF*8 B and big-loop geoms reach ~6 M DoF -> 300 vectors ~ 14 GB.
        # Paired with memoria_a_usar 44 G (job / 12 cores); measured peak ~31 G
        # at MaxSize 200, so ~40 G expected here -- inside 44 G, and SGE places
        # the job only where 44 G is free so a shared node can't oversubscribe.
        cfg['Solver']['Linear']['MaxSize'] = gmres_maxsize
        json.dump(cfg, open(mag._sim_config, 'w'), indent=2)

        n_tet = None
        try:
            import gmsh
            n_tet = len(gmsh.model.mesh.getElementsByType(4)[0])
        except Exception:
            pass

        flux_per_amp = float(np.array(mag.run()).ravel()[0])

        raw = PALACE_Inductance_Simulation.retrieve_magnetostatic_data_from_file(mag._output_data_dir)
        phi = float(np.ravel(raw['surface_Flux'])[0])
        i_inc = float(np.ravel(raw['terminal_I'])[0])

        n_ksp = converged = None
        logp = os.path.join(mag._output_data_dir, 'out.log')
        if os.path.exists(logp):
            txt = open(logp).read()
            n_ksp = sum('KSP residual' in ln for ln in txt.splitlines())
            converged = ('solver converged' in txt.lower()) and ('did not converge' not in txt.lower())

        row.update({
            'flux_per_amp': phi / i_inc,
            'Phi0_per_mA': (phi / i_inc) / PHI0 * 1e-3,
            'Phi_mag_Wb': phi, 'I_inc_A': i_inc,
            'n_tet': n_tet, 'n_ksp': n_ksp, 'converged': converged,
        })
    except Exception as e:
        row['error'] = f"{type(e).__name__}: {e}".replace('\n', ' ')[:300]
        traceback.print_exc()
        try:  # keep the Palace log of the failure for later diagnosis
            src = os.path.join(sim_dir, 'outputFiles', 'out.log')
            if os.path.exists(src):
                tag = '_'.join(f"{params[a]:g}" for a in AXIS_ORDER)
                fdir = os.path.join(os.path.dirname(workdir.rstrip('/')), 'failures')
                os.makedirs(fdir, exist_ok=True)
                shutil.copy(src, os.path.join(fdir, f"g{grow_um:g}_{tag}.log"))
        except Exception:
            pass
    finally:
        shutil.rmtree(sim_dir, ignore_errors=True)
        gc.collect()
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--shard', type=int, required=True)      # 1..nshards
    ap.add_argument('--nshards', type=int, required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--palace-bin', required=True)
    ap.add_argument('--num-cpus', type=int, default=8)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    _patch_gmsh_multipoly()

    pts = all_points()[args.shard - 1::args.nshards]
    csv_path = os.path.join(args.outdir, f'shard_{args.shard:03d}.csv')
    workdir = os.path.join(args.outdir, f'work_{args.shard:03d}')
    os.makedirs(workdir, exist_ok=True)

    done = set()
    if os.path.exists(csv_path):
        with open(csv_path) as f:
            for r in csv.DictReader(f):
                if r.get('flux_per_amp', '') == '':
                    continue          # error rows get retried on resume
                done.add(tuple(round(float(r[a]), 2) for a in AXIS_ORDER))
    new = not os.path.exists(csv_path)
    fh = open(csv_path, 'a', newline='')
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader(); fh.flush()

    print(f"shard {args.shard}/{args.nshards}: {len(pts)} points, {len(done)} already done", flush=True)
    for k, tup in enumerate(pts):
        params = {a: float(v) for a, v in zip(AXIS_ORDER, tup)}
        key = tuple(round(params[a], 2) for a in AXIS_ORDER)
        if key in done:
            continue
        row = run_one(params, args.palace_bin, args.num_cpus, workdir)
        w.writerow(row); fh.flush()
        print(f"  [{k+1}/{len(pts)}] {key} -> "
              f"{row['Phi0_per_mA'] if row['Phi0_per_mA'] != '' else row['error']}  "
              f"({row['seconds']}s)", flush=True)
    fh.close()
    shutil.rmtree(workdir, ignore_errors=True)
    print(f"shard {args.shard} DONE", flush=True)


if __name__ == '__main__':
    main()
