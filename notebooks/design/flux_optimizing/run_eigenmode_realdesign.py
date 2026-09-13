#!/usr/bin/env python
"""Eigenmode corroboration on the REAL target design (squid_array_coupling's
build_design.py: N_SQUIDS=23, squid_gap=25um, arm_width=4um, arm_semi=15um,
segment_length=6um, meandered readout at RES_LENGTH_MM=3.74575mm, their own
converged f_r=7.9355 GHz) -- as opposed to flux_optimizing's own 5-axis sweep
corners. Same method as run_eigenmode.py + extract_eigen_M.py (array metal IN,
loops left OPEN -- no junction ports, so L_J is irrelevant here and the only
new physics vs. a driven uniform-current model is the resonator's own
self-consistent current profile), combined into one script since this design's
geometry (fixed launch-point meander, not route_clearance-parameterised)
doesn't fit flux_optimizing's geometry_anchors() helper.

    python run_eigenmode_realdesign.py --outdir eigen/realdesign \
        --palace-bin /path/to/palace --num-cpus 32
"""
import argparse
import json
import os
import sys
import time
import traceback

os.environ.setdefault('QISKIT_METAL_HEADLESS', '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import shapely
from scipy.special import ellipk

SAC_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'squid_array_coupling')
sys.path.insert(0, os.path.abspath(SAC_DIR))
import build_design as sac_bd  # noqa: E402  -- squid_array_coupling's own module

PHI0 = 2.067833848e-15
MU0 = 4e-7 * np.pi
C0 = 2.99792458e8
EPS_R = 11.45
L0 = 1.0e-3
AREA_TO_SI = L0 ** 2
VOL_TO_SI = L0 ** 3
CPW_WIDTH_UM, CPW_GAP_UM = 20.0, 12.25          # sac_bd.build_design's cpw_width/cpw_gap
A_UM = CPW_WIDTH_UM / 2
B_UM = A_UM + CPW_GAP_UM
KMOD = A_UM / B_UM
KK = float(ellipk(KMOD ** 2))
SHORT_STANDOFF_UM = 15.0


def loop_aperture_polys_m(n_squids, squid_gap_um, arm_width_um, arm_semi_length_um,
                          josephson_gap_um, segment_length_um, array_pos_x_um,
                          array_pos_y_um, grow_um=0.3):
    """Same formula as flux_optimizing/build_design.loop_aperture_polys
    (copied inline -- the two build_design modules can't both be on sys.path
    under the same module name). Returns a shapely MultiPolygon, metres."""
    full_w = (arm_width_um + 2 * arm_semi_length_um + josephson_gap_um
             + arm_width_um + segment_length_um)
    loop_w = 2 * arm_semi_length_um + josephson_gap_um
    g = grow_um
    rects = []
    for i in range(int(n_squids)):
        x_left = array_pos_x_um + i * full_w + 0.5 * arm_width_um
        rects.append(shapely.geometry.box(
            (x_left - g) * 1e-6, (array_pos_y_um - 0.5 * squid_gap_um - g) * 1e-6,
            (x_left + loop_w + g) * 1e-6, (array_pos_y_um + 0.5 * squid_gap_um + g) * 1e-6))
    return shapely.MultiPolygon(rects)


def bz_cpw_to_I(bz_abs, b_um=B_UM):
    x_um = (A_UM + b_um) / 2
    sqrt_fac = np.sqrt((x_um ** 2 - A_UM ** 2) * (b_um ** 2 - x_um ** 2)) * 1e-12
    return bz_abs * 4 * KK * sqrt_fac / (MU0 * (b_um * 1e-6))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--palace-bin', required=True)
    ap.add_argument('--num-cpus', type=int, default=32)
    ap.add_argument('--number-of-freqs', type=int, default=6)
    ap.add_argument('--starting-freq', type=float, default=5e9)
    ap.add_argument('--solver-tol', type=float, default=1e-5)
    ap.add_argument('--solver-maxits', type=int, default=300)
    ap.add_argument('--fillet-resolution', type=int, default=12)
    ap.add_argument('--z-offset-um', type=float, default=3.0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation

    design = sac_bd.build_design()   # every default: N_SQUIDS=23, L_J=15nH
                                     # (irrelevant here -- no lumped ports),
                                     # RES_LENGTH_MM=3.74575, squid_gap=25um, ...
    design.rebuild()

    route = design.components['readout_1']
    route_len_um = route.length * 1e3
    v_ph = C0 / ((1 + EPS_R) / 2) ** 0.5
    f_analytic = v_ph / (4 * route_len_um * 1e-6)
    short = design.components['launch_point2']
    short_x_um = float(str(short.options.pos_x)[:-2])
    short_y_um = float(str(short.options.pos_y)[:-2])
    n_jj = len(design.qgeometry.tables['junction'].loc[
        design.qgeometry.tables['junction']['component']
        == design.components[sac_bd.ARRAY_COMPONENT_NAME].id])
    print(f"route length = {route_len_um:.0f} um (design RES_LENGTH_MM="
         f"{sac_bd.RES_LENGTH_MM}mm)   analytic bare lambda/4 f = "
         f"{f_analytic/1e9:.2f} GHz   (their own converged run: 7.9355 GHz)")
    print(f"short pin at ({short_x_um}, {short_y_um}) um   n_jj={n_jj} (left OPEN)",
         flush=True)

    apert = loop_aperture_polys_m(
        sac_bd.N_SQUIDS, sac_bd.SQUID_GAP_UM, sac_bd.ARM_WIDTH_UM,
        sac_bd.ARM_SEMI_LENGTH_UM, sac_bd.JOSEPHSON_GAP_UM, sac_bd.SEGMENT_LENGTH_UM,
        sac_bd.ARRAY_POS_X_UM, sac_bd.ARRAY_POS_Y_UM)
    n_loops = len(apert.geoms)
    apert_bounds_mm = [tuple(v * 1e3 for v in g.bounds) for g in apert.geoms]

    opts = {
        'mesh_refinement': 0, 'dielectric_material': 'silicon',
        'starting_freq': a.starting_freq, 'number_of_freqs': a.number_of_freqs,
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
        x1, y1, x2, y2 = apert.bounds
        pad = 8e-6
        eig.fine_mesh_in_rectangle(x1 - pad, y1 - pad, x2 + pad, y2 + pad,
                                   min_size=1.5e-6, max_size=60e-6, taper_dist_min=6e-6)
        eig.fine_mesh_components([sac_bd.ARRAY_COMPONENT_NAME], min_size=1.5e-6,
                                 max_size=60e-6, taper_dist_min=6e-6, metals_only=True)
        eig.fine_mesh_components(['readout_1'], min_size=6e-6, max_size=150e-6,
                                 taper_dist_min=40e-6, taper_dist_max=200e-6,
                                 metals_only=False)
        eig.prepare_simulation()
        eig.run()
    except Exception:
        traceback.print_exc()
        json.dump({'status': 'failed'}, open(os.path.join(a.outdir, 'result.json'), 'w'))
        sys.exit(1)

    OF = eig._output_data_dir
    eig_csv = os.path.join(OF, 'eig.csv')
    if not os.path.exists(eig_csv):
        print("FAILED: missing eig.csv")
        json.dump({'status': 'failed', 'reason': 'missing eig.csv'},
                  open(os.path.join(a.outdir, 'result.json'), 'w'))
        sys.exit(1)
    freqs = pd.read_csv(eig_csv).to_numpy()[:, 1]
    print(f"eigenfrequencies (GHz): {freqs}", flush=True)

    # ---- extraction (inline: this design's geometry doesn't fit
    # flux_optimizing's geometry_anchors(), so it can't reuse extract_eigen_M.py
    # as-is) ----
    def _read_pvtu(cyc, tries=4, delay_s=10):
        path = f"{OF}/paraview/eigenmode/Cycle{cyc:06d}/data.pvtu"
        last_err = None
        for attempt in range(tries):
            try:
                import pyvista as pv
                pv.OFF_SCREEN = True
                vol = pv.read(path).compute_cell_sizes(volume=True)
                vc = vol.point_data_to_cell_data()
                missing = [f for f in ('B_real', 'B_imag') if f not in vc.array_names]
                if missing:
                    raise KeyError(f"missing fields: {missing} (have {vc.array_names})")
                return vol, vc
            except Exception as e:
                last_err = e
                print(f"  [read attempt {attempt+1}/{tries}] {path}: {e} -- retry {delay_s}s",
                     flush=True)
                time.sleep(delay_s)
        raise RuntimeError(f"could not read {path} after {tries} attempts: {last_err}")

    def load_cycle(cyc):
        vol, vc = _read_pvtu(cyc)
        try:
            dV = np.asarray(vc['Volume']) * VOL_TO_SI
            Um_tot = float(np.sum(np.asarray(vc['U_m']) * dV))
        except Exception as e:
            print(f"  (U_m unavailable for cycle {cyc}: {e})")
            Um_tot = float('nan')
        sl = vol.slice(normal='z', origin=(0, 0, a.z_offset_um * 1e-3)).compute_cell_sizes(area=True)
        slc = sl.point_data_to_cell_data()
        sc = sl.cell_centers().points
        sBz = np.asarray(slc['B_real'])[:, 2] + 1j * np.asarray(slc['B_imag'])[:, 2]
        sA = np.asarray(sl['Area']) * AREA_TO_SI
        return sc, sBz, sA, Um_tot

    def total_loop_flux(sc, sBz, sA):
        tot = 0.0
        for (x0, y0, x1, y1) in apert_bounds_mm:
            m = (sc[:, 0] > x0) & (sc[:, 0] < x1) & (sc[:, 1] > y0) & (sc[:, 1] < y1)
            tot += abs(np.sum(sBz[m] * sA[m]))
        return tot

    n_cyc = len(freqs)
    scores, flux_scores = [], []
    for c in range(1, n_cyc + 1):
        try:
            sc, sBz, sA, _ = load_cycle(c)
            flux_scores.append(total_loop_flux(sc, sBz, sA))
            scores.append(abs(freqs[c - 1] - f_analytic / 1e9))
        except Exception as e:
            print(f"  cycle {c}: SKIPPED ({e})")
            flux_scores.append(-1.0)
            scores.append(np.inf)
    for c, (f, fl, s) in enumerate(zip(freqs, flux_scores, scores), 1):
        print(f"  cycle {c}: f={f:.3f} GHz  |f-f_analytic|={s:.3f} GHz   "
             f"total|Phi_loops|={fl:.3e} Wb")
    if min(scores) == np.inf:
        print("FAILED: no cycle usable")
        sys.exit(1)
    cyc = int(np.argmin(scores)) + 1
    if flux_scores[cyc - 1] < max(flux_scores) * 0.1:
        print(f"  WARNING: picked cycle {cyc} has much less array-flux-linkage "
             f"than another cycle")
    print(f"-> picked cycle {cyc} (f={freqs[cyc-1]:.4f} GHz)")

    f_r = freqs[cyc - 1] * 1e9
    sc, sBz, sA, Um_tot = load_cycle(cyc)
    try:
        dE = pd.read_csv(f"{OF}/domain-E.csv").to_numpy()
        emag = dE[cyc - 1, 2]
        print(f"cross-check: domain-E.csv E_mag={emag:.5e} J   U_m_tot/E_mag={Um_tot/emag:.4f}")
    except Exception as e:
        emag = np.nan
        print(f"cross-check skipped ({e})")

    Phi = np.zeros(n_loops, complex)
    for k, (x0, y0, x1, y1) in enumerate(apert_bounds_mm):
        m = (sc[:, 0] > x0) & (sc[:, 0] < x1) & (sc[:, 1] > y0) & (sc[:, 1] < y1)
        Phi[k] = complex(np.sum(sBz[m] * sA[m]))

    x_short_mm = (short_x_um - SHORT_STANDOFF_UM) / 1e3
    y_short_mm = short_y_um / 1e3

    def bz_at(x_mm, y_mm, r_mm=0.004):
        m = (np.abs(sc[:, 0] - x_mm) < r_mm) & (np.abs(sc[:, 1] - y_mm) < r_mm)
        return complex(np.mean(sBz[m])) if m.any() else np.nan

    vals = []
    for side in (-1, 1):
        yg = y_short_mm + side * ((A_UM + B_UM) / 2) * 1e-3
        bz = bz_at(x_short_mm, yg)
        if np.isfinite(bz):
            vals.append(bz_cpw_to_I(abs(bz)))
    if not vals:
        print("FAILED: could not sample B_z near the short")
        sys.exit(1)
    I_short = float(np.mean(vals))
    print(f"I_short (at ({x_short_mm*1e3:.0f},{y_short_mm*1e3:.0f})um, standoff "
         f"{SHORT_STANDOFF_UM}um) = {I_short:.4e} A ({len(vals)}/2 samples)")

    phi0g = np.angle(Phi[np.argmax(np.abs(Phi))])
    Phi_s = np.real(Phi * np.exp(-1j * phi0g))
    M_k = Phi_s / I_short
    M_total = float(np.sum(Phi_s) / I_short)

    print(f"\nf_r (cycle {cyc}) = {f_r/1e9:.4f} GHz")
    print(f"M_total (eigenmode, shorted-end ref) = {M_total*1e12:.3f} pH")

    result = {
        'status': 'ok', 'design': 'squid_array_coupling',
        'n_squids': sac_bd.N_SQUIDS, 'squid_gap_um': sac_bd.SQUID_GAP_UM,
        'arm_width_um': sac_bd.ARM_WIDTH_UM, 'arm_semi_length_um': sac_bd.ARM_SEMI_LENGTH_UM,
        'segment_length_um': sac_bd.SEGMENT_LENGTH_UM, 'res_length_mm': sac_bd.RES_LENGTH_MM,
        'route_len_um': route_len_um, 'f_analytic_hz': f_analytic,
        'freqs_ghz': freqs.tolist(), 'cycle': cyc, 'f_r_hz': f_r,
        'M_total_H': M_total, 'I_short_A': I_short,
        'Um_tot_J': Um_tot, 'E_mag_domainE_J': float(emag),
        'seconds': round(time.time() - t0, 1),
    }
    json.dump(result, open(os.path.join(a.outdir, 'result.json'), 'w'), indent=2)
    np.savez(os.path.join(a.outdir, f'eigen_M_cyc{cyc}.npz'),
            loop=np.arange(n_loops), Phi_signed=Phi_s, M_k=M_k, **result)

    import thesis_style as ts
    ts.use()
    fig, ax = plt.subplots(2, 1, figsize=ts.figsize(1.5, ratio=0.6), sharex=True)
    ax[0].bar(range(n_loops), Phi_s * 1e15, color=ts.COLORS['blue'])
    ax[0].set_title(f'$\\Phi_k$ (fWb)   z={a.z_offset_um:.1f}um above plane')
    ax[1].bar(range(n_loops), M_k * 1e12, color=ts.COLORS['teal'])
    ax[1].set_title('$M_k = \\Phi_k / I_{short}$  (pH)')
    ax[1].set_xlabel('loop index')
    fig.suptitle(f"squid_array_coupling real design, cycle {cyc}, $f_r$={f_r/1e9:.3f} GHz\n"
                f"$M_{{total}}$={M_total*1e12:.2f} pH", fontsize=9)
    ts.finalize(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, f'eigen_M_cyc{cyc}.png'), dpi=160, facecolor='white')
    print(f"done ({time.time()-t0:.0f}s) -> {a.outdir}/result.json")


if __name__ == '__main__':
    main()
