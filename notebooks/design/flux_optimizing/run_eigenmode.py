#!/usr/bin/env python
"""Eigenmode corroboration for ONE point of the flux-threading sweep.

Same design as sweep_flux.py (array metal IN, loops left OPEN -- no junction
ports, exactly matching the magnetostatic model so the only new physics is the
resonator's own self-consistent standing-wave current, not loop screening),
but solved as a Palace EIGENMODE problem instead of driven-magnetostatic: the
readout hairpin's own quarter-wave mode threads the array, and we read the
per-loop flux directly off that mode's B-field (see extract_eigen_M.py).

Unlike the magnetostatic sweep, this uses the REAL substrate thickness
(build_design's default CHIP_SIZE_Z_UM=500um, not the 200um the magnetostatic
sweep used to save mesh volume) since it sets the real f_r here.

    python run_eigenmode.py --route-clearance-um 30 --segment-length-um 3 \
        --squid-gap-um 40 --arm-semi-length-um 24 --arm-width-um 8 \
        --outdir eigen/best --palace-bin /path/to/palace --num-cpus 16
"""
import argparse
import json
import os
import sys
import time
import traceback

os.environ.setdefault('QISKIT_METAL_HEADLESS', '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

EPS_R = 11.45          # silicon, thesis convention (transmon_coupler/*.py)
C0 = 2.99792458e8      # m/s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--route-clearance-um', type=float, required=True)
    ap.add_argument('--segment-length-um', type=float, required=True)
    ap.add_argument('--squid-gap-um', type=float, required=True)
    ap.add_argument('--arm-semi-length-um', type=float, required=True)
    ap.add_argument('--arm-width-um', type=float, required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--palace-bin', required=True)
    ap.add_argument('--num-cpus', type=int, default=16)
    ap.add_argument('--number-of-freqs', type=int, default=4)
    ap.add_argument('--starting-freq', type=float, default=None,
                    help='Hz. Default: 0.6 * analytic bare lambda/4 estimate '
                         '(v_ph = c/sqrt((1+eps_Si)/2), f = v_ph/(4*route_length)).')
    ap.add_argument('--solver-tol', type=float, default=1e-5)
    ap.add_argument('--solver-maxits', type=int, default=300)
    ap.add_argument('--fillet-resolution', type=int, default=12)
    ap.add_argument('--array-ground-margin-um', type=float, default=None,
                    help='None = build_design default (6um rectangular pocket).')
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    from build_design import build_design, ARRAY_COMPONENT_NAME
    from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation

    params = {k: getattr(a, k) for k in
             ('route_clearance_um', 'segment_length_um', 'squid_gap_um',
              'arm_semi_length_um', 'arm_width_um')}
    print(f"config: {params}", flush=True)

    design = build_design(
        chip_margin_um=250, include_array_metal=True, ground_near_end=True,
        array_ground_margin_um=a.array_ground_margin_um,   # None -> the new 6um rectangular default
        route_clearance_um=a.route_clearance_um,
        segment_length_um=a.segment_length_um,
        squid_gap_um=a.squid_gap_um,
        arm_semi_length_um=a.arm_semi_length_um,
        arm_width_um=a.arm_width_um)
    # NOTE: chip_size_z_um NOT overridden -- keep build_design's real 500um
    # substrate default (the magnetostatic sweep shrank it to 200um purely to
    # cut mesh volume; here it sets the real f_r).
    design.rebuild()

    route = design.components['readout_1']
    route_len_um = route.length * 1e3
    v_ph = C0 / ((1 + EPS_R) / 2) ** 0.5
    f_analytic = v_ph / (4 * route_len_um * 1e-6)
    starting_freq = a.starting_freq if a.starting_freq is not None else 0.6 * f_analytic
    print(f"route length = {route_len_um:.0f} um   analytic bare lambda/4 f "
         f"= {f_analytic/1e9:.2f} GHz   starting_freq = {starting_freq/1e9:.2f} GHz",
         flush=True)

    n_jj = len(design.qgeometry.tables['junction'].loc[
        design.qgeometry.tables['junction']['component']
        == design.components[ARRAY_COMPONENT_NAME].id])
    print(f"array junction rows: {n_jj} (left OPEN -- no lumped ports, "
         f"matching the magnetostatic sweep)", flush=True)

    opts = {
        'mesh_refinement': 0, 'dielectric_material': 'silicon',
        'starting_freq': starting_freq, 'number_of_freqs': a.number_of_freqs,
        'solns_to_save': a.number_of_freqs,
        'solver_order': 1, 'solver_tol': a.solver_tol,
        'solver_maxits': a.solver_maxits, 'fillet_resolution': a.fillet_resolution,
        'palace_dir': a.palace_bin, 'num_cpus': a.num_cpus,
    }
    t0 = time.time()
    try:
        eig = PALACE_Eigenmode_Simulation(
            name='eig', metal_design=design, sim_parent_directory=a.outdir,
            mode='simPC', meshing='GMSH', user_options=opts, create_files=True)
        eig.add_metallic(1)
        eig.add_ground_plane()

        # Fine mesh over the loop-aperture band (same recipe as sweep_flux.py's
        # magnetostatic solve -- 1.5um resolves the flux-integration cuts) PLUS
        # the resonator trace/CPW-gap, which now also carries real mode physics.
        from build_design import loop_aperture_polys
        int_area = loop_aperture_polys(
            squid_gap_um=a.squid_gap_um, arm_semi_length_um=a.arm_semi_length_um,
            arm_width_um=a.arm_width_um, segment_length_um=a.segment_length_um)
        x1, y1, x2, y2 = int_area.bounds
        pad = 8e-6
        eig.fine_mesh_in_rectangle(x1 - pad, y1 - pad, x2 + pad, y2 + pad,
                                   min_size=1.5e-6, max_size=60e-6, taper_dist_min=6e-6)
        eig.fine_mesh_components([ARRAY_COMPONENT_NAME], min_size=1.5e-6, max_size=60e-6,
                                 taper_dist_min=6e-6, metals_only=True)
        eig.fine_mesh_components(['readout_1'], min_size=6e-6, max_size=150e-6,
                                 taper_dist_min=40e-6, taper_dist_max=200e-6,
                                 metals_only=False)
        eig.prepare_simulation()
        eig.run()
    except Exception:
        traceback.print_exc()
        json.dump({'status': 'failed', **params, 'route_len_um': route_len_um},
                  open(os.path.join(a.outdir, 'result.json'), 'w'), indent=2)
        sys.exit(1)

    eig_csv = os.path.join(eig._output_data_dir, 'eig.csv')
    if not os.path.exists(eig_csv):
        print(f"FAILED: run() returned but {eig_csv} missing (killed mid-run?)")
        json.dump({'status': 'failed', 'reason': 'missing eig.csv', **params,
                  'route_len_um': route_len_um},
                  open(os.path.join(a.outdir, 'result.json'), 'w'), indent=2)
        sys.exit(1)

    import pandas as pd
    freqs = pd.read_csv(eig_csv).to_numpy()[:, 1]
    print(f"eigenfrequencies (GHz): {freqs}", flush=True)
    json.dump({
        'status': 'ok', **params, 'route_len_um': route_len_um,
        'f_analytic_hz': f_analytic, 'starting_freq_hz': starting_freq,
        'freqs_ghz': freqs.tolist(), 'output_data_dir': eig._output_data_dir,
        'seconds': round(time.time() - t0, 1),
    }, open(os.path.join(a.outdir, 'result.json'), 'w'), indent=2)
    print(f"done ({time.time()-t0:.0f}s) -> {a.outdir}/result.json", flush=True)


if __name__ == '__main__':
    main()
