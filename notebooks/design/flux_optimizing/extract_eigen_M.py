#!/usr/bin/env python
"""Per-loop mutual inductance from an eigenmode B-field, corroborating the
magnetostatic flux sweep.

Method (same field-extraction technique as
../squid_array_coupling/extract_M_bfield.py, no junction ports here -- the
array loops are left OPEN, exactly matching the magnetostatic sweep's model,
so the only new physics vs. that sweep is the resonator's own self-consistent
current profile, not loop screening):

  1. Flux through loop k: integral of B_z over the loop's aperture (the SAME
     ``loop_aperture_polys`` box the magnetostatic sweep integrated), on a
     horizontal cut a few um ABOVE the metal plane (on z=0 B.n=0 by the PEC BC).

  2. Reference current: the exact zero-thickness-PEC CPW inversion (Lahl &
     Woerdenweber, APL 81 505 (2002); Bothner et al. arXiv:1202.5963),
        B_z,gap(x) = mu0 I b / (4 K(k) sqrt((x^2-a^2)(b^2-x^2))),  a=w/2, b=a+gap
     sampled in the CPW gap a short standoff BACK from the ShortToGround pin
     (the resonator's current antinode) -- NOT the local current at each loop.
     This is the corroboration-specific choice (user, 2026-09-04): every loop's
     M is referenced to the SAME point (the short), so the M_k are directly
     comparable to each other and to M = Sum Phi_k / I_probe from the
     magnetostatic sweep (which used a uniform injected current, i.e. every
     point on the line implicitly "is" the reference).

  3. M_k = Phi_k / I_short (H).  M_total = Sum_k Phi_k / I_short, directly
     comparable to sweep_flux.py's Phi0_per_mA * PHI0 * 1e-3.

UNITS: Palace mesh L0 = 1e-3 (mesh coords in mm); B (T) and U_m (J/m^3) in SI.
    area_SI [m^2] = area_mesh [mm^2] * 1e-6

Usage:  python extract_eigen_M.py <run_outdir> [--cycle N] [--z-offset-um 3.0]
    <run_outdir> is the --outdir passed to run_eigenmode.py (holds result.json).
    --cycle: which of the solved eigenmodes to analyze. Default: auto-picks the
    mode with the largest total |flux| linked to the array (the resonator
    fundamental, as opposed to a spurious box/substrate mode).
Writes: <run_outdir>/eigen_M_cyc<C>.npz, <run_outdir>/eigen_M_cyc<C>.png
"""
import argparse
import json
import os
import sys

os.environ.setdefault('QISKIT_METAL_HEADLESS', '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import numpy as np
import pandas as pd
import pyvista as pv
from scipy.special import ellipk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import geometry_anchors, loop_aperture_polys  # noqa: E402

pv.OFF_SCREEN = True

MU0 = 4e-7 * np.pi
PHI0 = 2.067833848e-15
L0 = 1.0e-3
AREA_TO_SI = L0 ** 2
VOL_TO_SI = L0 ** 3

CPW_WIDTH_UM = 20.0    # build_design defaults -- never varied in this sweep
CPW_GAP_UM = 12.25
A_UM = CPW_WIDTH_UM / 2
B_UM = A_UM + CPW_GAP_UM
KMOD = A_UM / B_UM
KK = float(ellipk(KMOD ** 2))
SHORT_STANDOFF_UM = 15.0   # sample this far back from the short, along the
                          # straight lead (RouteAnchors lead end_straight=30um)


def bz_cpw_to_I(bz_abs, b_um=B_UM):
    """Invert the exact zero-thickness-PEC CPW gap field -> current (A)."""
    x_um = (A_UM + b_um) / 2
    sqrt_fac = np.sqrt((x_um ** 2 - A_UM ** 2) * (b_um ** 2 - x_um ** 2)) * 1e-12
    return bz_abs * 4 * KK * sqrt_fac / (MU0 * (b_um * 1e-6))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_outdir')
    ap.add_argument('--cycle', type=int, default=None)
    ap.add_argument('--z-offset-um', type=float, default=3.0)
    a = ap.parse_args()

    result = json.load(open(os.path.join(a.run_outdir, 'result.json')))
    if result.get('status') != 'ok':
        print(f"run did not succeed: {result}")
        sys.exit(1)
    params = {k: result[k] for k in ('route_clearance_um', 'segment_length_um',
                                    'squid_gap_um', 'arm_semi_length_um', 'arm_width_um')}
    OF = result['output_data_dir']
    print(f"config: {params}")
    print(f"f_analytic = {result['f_analytic_hz']/1e9:.2f} GHz   "
         f"solved freqs (GHz) = {result['freqs_ghz']}")

    ga = geometry_anchors(route_clearance_um=params['route_clearance_um'],
                          segment_length_um=params['segment_length_um'],
                          squid_gap_um=params['squid_gap_um'],
                          arm_semi_length_um=params['arm_semi_length_um'],
                          arm_width_um=params['arm_width_um'])
    apert = loop_aperture_polys(squid_gap_um=params['squid_gap_um'],
                                arm_semi_length_um=params['arm_semi_length_um'],
                                arm_width_um=params['arm_width_um'],
                                segment_length_um=params['segment_length_um'])
    n_loops = len(apert.geoms)
    apert_bounds_mm = [tuple(v * 1e3 for v in g.bounds) for g in apert.geoms]  # m -> mm

    freqs = np.asarray(result['freqs_ghz'])
    n_cyc = len(freqs)

    def _read_pvtu(cyc, tries=4, delay_s=10):
        # Parallel-piece VTU reads have occasionally come back with a field
        # missing (e.g. U_m) when hit shortly after the writing MPI job
        # finished, on a compute node other than the one that wrote it --
        # looks like NFS attribute-cache staleness across nodes (a plain
        # re-read from the frontend right after has always had every field).
        # Retry with backoff rather than fail the whole extraction over it.
        path = f"{OF}/paraview/eigenmode/Cycle{cyc:06d}/data.pvtu"
        last_err = None
        for attempt in range(tries):
            try:
                vol = pv.read(path).compute_cell_sizes(volume=True)
                vc = vol.point_data_to_cell_data()
                missing = [f for f in ('B_real', 'B_imag') if f not in vc.array_names]
                if missing:
                    raise KeyError(f"missing fields after read: {missing} "
                                  f"(have: {vc.array_names})")
                return vol, vc
            except Exception as e:
                last_err = e
                print(f"  [read attempt {attempt+1}/{tries}] {path}: {e} "
                     f"-- retrying in {delay_s}s", flush=True)
                import time
                time.sleep(delay_s)
        raise RuntimeError(f"could not read a complete {path} after {tries} "
                           f"attempts: {last_err}")

    def load_cycle(cyc):
        vol, vc = _read_pvtu(cyc)
        try:
            dV = np.asarray(vc['Volume']) * VOL_TO_SI
            Um_tot = float(np.sum(np.asarray(vc['U_m']) * dV))
        except Exception as e:
            print(f"  (U_m cross-check unavailable for cycle {cyc}: {e} -- "
                 f"not required for M_k, continuing without it)")
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

    if a.cycle is not None:
        cyc = a.cycle
    else:
        # auto-pick the mode closest to the analytic bare-lambda/4 estimate --
        # print all candidates (frequency AND array-flux-linkage) so the
        # choice is auditable. Palace doesn't always write every requested
        # mode to paraview; skip cycles whose fields are genuinely absent
        # rather than fail the whole run over a mode we don't need.
        #
        # NOT total|Phi_loops| (tried first, wrong for the big-loop "best"
        # config): a large SQUID array is itself a resonant structure (loop
        # apertures + nearby ground), and one of ITS own modes can couple to
        # the flux-integration boxes far more strongly than the actual
        # readout-line lambda/4 mode -- e.g. for clr30/seg3/gap40/semi24/aw8,
        # cycle 1 (8.742 GHz, matching f_analytic=8.74 GHz almost exactly)
        # had total|Phi_loops|=6.5e-12 Wb, but cycle 4 (29.35 GHz, a spurious
        # array/box mode) had 4.1e-11 Wb -- 6x more -- and would have been
        # picked, and reported a meaningless M. Frequency proximity to the
        # independently-computed bare-line estimate is a far more reliable
        # resonator-identification criterion than array-coupling strength.
        f_analytic_ghz = result['f_analytic_hz'] / 1e9
        scores, flux_scores = [], []
        for c in range(1, n_cyc + 1):
            try:
                sc, sBz, sA, _ = load_cycle(c)
                flux_scores.append(total_loop_flux(sc, sBz, sA))
                scores.append(abs(freqs[c - 1] - f_analytic_ghz))
            except Exception as e:
                print(f"  cycle {c}: SKIPPED (no usable field data: {e})")
                flux_scores.append(-1.0)
                scores.append(np.inf)
        for c, (f, fl, s) in enumerate(zip(freqs, flux_scores, scores), 1):
            print(f"  cycle {c}: f={f:.3f} GHz  |f-f_analytic|={s:.3f} GHz   "
                 f"total|Phi_loops|={fl:.3e} Wb (mode-normalised)")
        if min(scores) == np.inf:
            print("FAILED: no cycle had usable field data")
            sys.exit(1)
        cyc = int(np.argmin(scores)) + 1
        if flux_scores[cyc - 1] < max(flux_scores) * 0.1:
            print(f"  WARNING: picked cycle {cyc} has much less array-flux-"
                 f"linkage than another cycle -- double check with --cycle")
        print(f"-> auto-picked cycle {cyc} (f={freqs[cyc-1]:.3f} GHz, "
             f"closest to f_analytic={f_analytic_ghz:.2f} GHz)")

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

    # ---- reference current: sampled at the SHORTED end, not per loop ----
    x_short_mm = (ga['term_x'] - SHORT_STANDOFF_UM) / 1e3
    y_bot_mm = ga['pass_y_bot'] / 1e3

    def bz_at(x_mm, y_mm, r_mm=0.004):
        m = (np.abs(sc[:, 0] - x_mm) < r_mm) & (np.abs(sc[:, 1] - y_mm) < r_mm)
        return complex(np.mean(sBz[m])) if m.any() else np.nan

    vals = []
    for side in (-1, 1):
        yg = y_bot_mm + side * ((A_UM + B_UM) / 2) * 1e-3
        bz = bz_at(x_short_mm, yg)
        if np.isfinite(bz):
            vals.append(bz_cpw_to_I(abs(bz)))
    if not vals:
        print("FAILED: could not sample B_z near the short -- check mesh/geometry")
        sys.exit(1)
    I_short = float(np.mean(vals))
    print(f"I_short (at shorted end, standoff {SHORT_STANDOFF_UM}um) = {I_short:.4e} A "
         f"(mode-normalised, {len(vals)}/2 gap-side samples)")

    phi0g = np.angle(Phi[np.argmax(np.abs(Phi))])
    Phi_s = np.real(Phi * np.exp(-1j * phi0g))
    M_k = Phi_s / I_short                     # H, signed, per loop
    M_total = float(np.sum(Phi_s) / I_short)  # H

    # magnetostatic comparison, if the sweep result for this exact config is on disk
    M_magneto = None
    for cand in ('sweep/sweep_results.csv', '../sweep/sweep_results.csv'):
        if os.path.exists(cand):
            df = pd.read_csv(cand)
            m = np.ones(len(df), bool)
            for k, v in params.items():
                m &= np.isclose(df[k], v)
            if m.any():
                M_magneto = float(df.loc[m, 'Phi0_per_mA'].iloc[0]) * 1e-3 * PHI0
            break

    print(f"\nf_r (cycle {cyc}) = {f_r/1e9:.4f} GHz")
    print(f"M_total (eigenmode, shorted-end ref) = {M_total*1e12:.3f} pH")
    if M_magneto is not None:
        print(f"M (magnetostatic, uniform-I ref)     = {M_magneto*1e12:.3f} pH   "
             f"ratio eigen/magneto = {M_total/M_magneto:.3f}")
    else:
        print("(no matching magnetostatic sweep row found on disk for comparison)")

    npz = os.path.join(a.run_outdir, f'eigen_M_cyc{cyc}.npz')
    np.savez(npz, loop=np.arange(n_loops), Phi_complex=Phi, Phi_signed=Phi_s,
            M_k=M_k, I_short_A=I_short, M_total_H=M_total,
            M_magneto_H=M_magneto if M_magneto is not None else np.nan,
            f_r_Hz=f_r, Um_tot_J=Um_tot, E_mag_domainE_J=emag,
            cycle=cyc, z_offset_um=a.z_offset_um, **params)
    print(f"saved {npz}")

    import thesis_style as ts
    ts.use()
    fig, ax = plt.subplots(2, 1, figsize=ts.figsize(1.5, ratio=0.6), sharex=True)
    ax[0].bar(range(n_loops), Phi_s * 1e15, color=ts.COLORS['blue'])
    ax[0].set_title(f'$\\Phi_k$ (fWb)   z={a.z_offset_um:.1f}um above plane')
    ax[1].bar(range(n_loops), M_k * 1e12, color=ts.COLORS['teal'])
    ax[1].set_title('$M_k = \\Phi_k / I_{short}$  (pH)')
    ax[1].set_xlabel('loop index')
    ttl = (f"{a.run_outdir}  cycle {cyc}  $f_r$={f_r/1e9:.3f} GHz\n"
          f"$M_{{total}}$={M_total*1e12:.2f} pH")
    if M_magneto is not None:
        ttl += f"   $M_{{magnetostatic}}$={M_magneto*1e12:.2f} pH  ratio={M_total/M_magneto:.3f}"
    fig.suptitle(ttl, fontsize=9)
    ts.finalize(ax)
    fig.tight_layout()
    png = os.path.join(a.run_outdir, f'eigen_M_cyc{cyc}.png')
    fig.savefig(png, dpi=160, facecolor='white')
    print(f"saved {png}")


if __name__ == '__main__':
    main()
