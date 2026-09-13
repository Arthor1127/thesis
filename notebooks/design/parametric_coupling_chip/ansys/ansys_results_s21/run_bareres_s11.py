"""ANSYS HFSS driven-S11 cross-check of PALACE's eigenmode-Qc method -- see
../HANDOFF.md.

Geometry (from ../geometry/bare_resonator.gds): a straight shorted quarter-wave
CPW resonator, w = 20 um, gap = 11.5 um, length ~3.773 mm on 500 um silicon
(er = 11.45).  One end is a real PEC short; the open end carries a single
**lumped port of resistance R**.  Lossless dielectric + PEC metal, so the only
loss channel is the port -> the loaded Q measured from the driven S11 IS Q_c.

Q_c is read two ways from the same sweep:
  * group delay of S11 at resonance:  Q_c = pi * f_r * tau_g(f_r) / 2
  * a complex one-port (Probst-style) resonance fit of S11(f)

Run (in the qiskit_qc env):
  python run_bareres_s11.py            # all three R values
  python run_bareres_s11.py 1e6        # just one

Env: HFSS_CORES (default 6), AEDT_VERSION (2024.2)
"""
import argparse
import json
import os
import sys
import time
import traceback

import numpy as np

os.environ.setdefault("SI_TANDELTA", "1.2e-5")   # small known dielectric loss so
#   the driven adaptive mesh converges on the resonance; Q_i is removed in post

HERE = os.path.dirname(os.path.abspath(__file__))
SHARED = os.path.abspath(os.path.join(HERE, "..", "..", "ansys_impedance_handoff",
                                      "ansys_results"))
GDS = os.path.abspath(os.path.join(HERE, "..", "geometry", "bare_resonator.gds"))
sys.path.insert(0, SHARED)
from gds_geometry import read_metal_polygons, bbox  # noqa: E402
import run_ansys_tests as R  # shared HFSS helpers  # noqa: E402

AEDT_VERSION = os.environ.get("AEDT_VERSION", "2024.2")
CORES = int(os.environ.get("HFSS_CORES", "6"))

ER_SILICON = 11.45
SUB_THICK_UM = 500.0
CPW_W_UM = 20.0
CPW_G_UM = 11.5
Z0_ANALYTIC = 49.961

R_VALUES = [2e5, 1e6, 5e6]
ANALYTIC_QC = {2e5: 3144, 1e6: 15720, 5e6: 78601}
PALACE_QC = {2e5: 3841.1, 1e6: 18495.4, 5e6: 78181.5}


# ---------------------------------------------------------------------------

def _rect(m, x0, x1, y0, y1, name):
    return m.create_polyline(points=[[x0, y0, 0], [x1, y0, 0], [x1, y1, 0], [x0, y1, 0]],
                             name=name, cover_surface=True, close_surface=True).name


def build(h, y_half_um=200.0, h_vac_um=200.0, port_len_um=10.0):
    """Straight CPW: PEC short at +x, lumped-port gap(s) at -x.

    The lumped element (driven port or eigenmode RLC) bridges a SINGLE CPW gap
    -- trace to ground -- so the resistance seen by the CPW mode is exactly the
    nominal R (a sheet spanning both gaps puts two R's in parallel = R/2).
    """
    m = h.modeler
    m.model_units = "um"
    p = read_metal_polygons(GDS)
    tb = bbox(p["trace"])                     # trace bbox (um)
    x_open, x_short = tb[0], tb[2]            # open end = min x, short end = max x
    w = tb[3] - tb[1]
    yt = w / 2.0
    yg = w / 2.0 + CPW_G_UM
    Y = y_half_um
    Lp = port_len_um

    sil = "silicon_er1145"
    if sil not in h.materials.material_keys:
        mat = h.materials.add_material(sil)
        mat.permittivity = ER_SILICON
        mat.dielectric_loss_tangent = float(os.environ.get("SI_TANDELTA", "0.0"))

    sub = m.create_box([x_open - 50, -Y, -SUB_THICK_UM],
                       [(x_short + 50) - (x_open - 50), 2 * Y, SUB_THICK_UM],
                       name="substrate", material=sil)
    m.create_box([x_open - 50, -Y, 0.0],
                 [(x_short + 50) - (x_open - 50), 2 * Y, h_vac_um],
                 name="vacuum", material="vacuum")

    # trace runs the full physical length, open at -x; PEC bridge shorts +x
    trace = _rect(m, x_open, x_short - 20, -yt, yt, "trace")
    gnd_l = _rect(m, x_open - 50, x_short, -Y, -yg, "gnd_l")
    gnd_r = _rect(m, x_open - 50, x_short, yg, Y, "gnd_r")
    short = _rect(m, x_short - 20, x_short, -yg, yg, "short_bridge")
    gnd = m.unite(["gnd_l", "gnd_r", "short_bridge"])
    gnd = gnd if isinstance(gnd, str) else "gnd_l"

    h.assign_perfecte_to_sheets(["trace"], name="PerfE_signal")
    h.assign_perfecte_to_sheets([gnd], name="PerfE_gnd")
    h.assign_perfecte_to_sheets([sub.bottom_face_z.id], name="PerfE_backplane")

    # lumped-element sheet: last Lp of the trace's open end, ONE gap (y in
    # [yt, yg]).  Touches the trace top edge and the ground inner edge.
    ps = _rect(m, x_open, x_open + Lp, yt, yg, "portsheet")
    x_tr0 = x_open

    # (a blanket length-mesh on these long thin CPW sheets hangs the mesher on
    # this machine -- disabled.  Instead the adaptive solve is forced through
    # many passes and driven off-resonance, and tan_delta makes dS non-blind.)
    ml = os.environ.get("BR_MESH", "")
    if ml:
        if not ml.endswith(("um", "mm")):
            ml = ml + "um"
        try:
            h.mesh.assign_length_mesh(["trace", gnd], inside_selection=False,
                                      maximum_length=ml, maximum_elements=600000,
                                      name="seed_cpw")
        except Exception as e:  # noqa: BLE001
            print("   mesh seed warning:", repr(e), flush=True)

    return dict(port_sheet=ps, x_open=x_open, x_tr0=x_tr0, yt=yt, yg=yg,
                x_short=x_short, x_line=x_open + Lp / 2.0, fr_guess_ghz=8.0)


def _iline(geo):
    """integration line across the one gap: trace edge -> ground edge."""
    x = geo["x_line"]
    return [[x, geo["yt"], 0.0], [x, geo["yg"], 0.0]]


def add_lumped_port(h, geo, Rohm):
    return h.lumped_port(assignment=geo["port_sheet"], integration_line=_iline(geo),
                         impedance=float(Rohm), renormalize=False, deembed=False,
                         name="p1")


# ---------------------------------------------------------------------------
# Q_c extraction
# ---------------------------------------------------------------------------

def _group_delay_qc(f, s11):
    ph = np.unwrap(np.angle(s11))
    tau_g = -np.gradient(ph, f) / (2 * np.pi)
    j = int(np.argmax(np.abs(tau_g - np.median(tau_g))))
    fr = f[j]
    # local parabolic peak of |tau_g|
    lo, hi = max(0, j - 8), min(len(f), j + 9)
    try:
        c = np.polyfit(f[lo:hi], tau_g[lo:hi], 2)
        fr = -c[1] / (2 * c[0])
        tau_peak = np.polyval(c, fr)
    except Exception:
        tau_peak = tau_g[j]
    return dict(fr_hz=float(fr), tau_g_peak_s=float(tau_peak),
                Qc=float(np.pi * fr * abs(tau_peak) / 2.0))


def _oneport_fit(f, s11):
    """Complex one-port (reflection) resonance fit -- Probst et al. 2015 style:

        S11 = env(f) * [ 1  -  (2 Ql / |Qc| e^{-i phi}) / (1 + 2 i Ql (f/fr - 1)) ]

    env(f) = a e^{i alpha} e^{-2 pi i (f - f0) tau}  is the cable/background.
    Returns Ql (loaded), Qc (external/coupling, from |Qc| and phi) and
    Qi = 1/(1/Ql - Re(1/Qc_complex)).
    """
    from scipy.optimize import least_squares
    f0 = f.mean()
    mag = np.abs(s11)
    ph = np.unwrap(np.angle(s11))
    # seed fr and Ql from the group-delay peak (robust for a 2 pi phase wrap)
    tau = -np.gradient(ph, f) / (2 * np.pi)
    j = int(np.argmax(np.abs(tau - np.median(tau))))
    fr0 = float(f[j])
    # FWHM of the |tau| peak -> loaded Q
    pk = abs(tau[j] - np.median(tau))
    inw = f[np.abs(tau - np.median(tau)) > 0.5 * pk]
    bw = (inw.max() - inw.min()) if inw.size > 2 else np.ptp(f) / 50
    ql0 = float(np.clip(fr0 / max(bw, f[1] - f[0]), 30, 5e6))

    def model(pp):
        fr, ql, qc_abs, phi, a, alpha, tau = pp
        env = a * np.exp(1j * alpha) * np.exp(-2j * np.pi * (f - f0) * tau)
        return env * (1 - (2 * ql / qc_abs) * np.exp(-1j * phi)
                      / (1 + 2j * ql * (f / fr - 1)))

    def resid(pp):
        d = model(pp) - s11
        return np.concatenate([d.real, d.imag])

    best = None
    for qcm in (1.0, 1.3, 3.0, 10.0):
        p0 = [fr0, ql0, ql0 * qcm, 0.0, np.median(mag),
              float(np.angle(s11[np.argmin(np.abs(f - f0))])), 0.0]
        try:
            sol = least_squares(resid, p0, method="trf", x_scale="jac",
                                max_nfev=60000)
            if best is None or sol.cost < best.cost:
                best = sol
        except Exception:
            pass
    if best is None:
        return {"error": "fit failed"}
    fr, ql, qc_abs, phi, a, alpha, tau = best.x
    ql, qc_abs = abs(ql), abs(qc_abs)
    qc_complex = qc_abs * np.exp(1j * phi)
    re_inv_qc = np.real(1.0 / qc_complex)
    qi = 1.0 / (1.0 / ql - re_inv_qc) if (1.0 / ql - re_inv_qc) > 0 else None
    return dict(fr_hz=float(fr), Ql=float(ql), Qc=float(qc_abs),
                phi_rad=float(phi), Qi=(float(qi) if qi else None),
                fit_cost=float(best.cost))


# ---------------------------------------------------------------------------

def run_one(Rohm, outdir, graphical, f0_ghz, f1_ghz, npts, fs_ghz, kind):
    tag = f"bareres_R{Rohm:.0e}".replace("+", "")
    print(f"\n{'='*60}\n  R = {Rohm:.3g} Ohm   analytic Qc = {ANALYTIC_QC.get(Rohm,'?')}"
          f"   (PALACE eig {PALACE_QC.get(Rohm,'?')})\n{'='*60}", flush=True)
    h = R.open_hfss(f"ansys_{tag}", graphical)
    try:
        geo = build(h)
        add_lumped_port(h, geo, Rohm)
        R.check_validation(h)
        h.save_project()

        # Adapt OFF the resonance -- for a uniform CPW line the transverse mode
        # is frequency-independent, so an off-resonance solve still resolves the
        # cross-section (+ the mesh seed), while dodging the near-singular FEM
        # system you get by solving exactly at a high-Q resonance.
        f_adapt = float(os.environ.get("BR_FADAPT", "8.30"))
        setup = h.create_setup("Setup1")
        setup.props["Frequency"] = f"{f_adapt}GHz"
        setup.props["MaximumPasses"] = int(os.environ.get("BR_MAXPASS", "25"))
        setup.props["MinimumPasses"] = int(os.environ.get("BR_MINPASS", "16"))
        setup.props["MaxDeltaS"] = float(os.environ.get("BR_DS", "1e-4"))
        setup.props["MinimumConvergedPasses"] = 2
        setup.update()
        # narrow window around the resonance; an interpolating (rational-model)
        # sweep solves at ~10-30 points (none on the singularity) and captures
        # the pole
        qexp = ANALYTIC_QC.get(Rohm) or (np.pi / 4) * (Rohm / Z0_ANALYTIC)
        lw = fs_ghz / qexp                       # linewidth in GHz
        span = float(os.environ.get("BR_SPAN_LW", "60"))
        f0_ghz = fs_ghz - span * lw
        f1_ghz = fs_ghz + span * lw
        kind = os.environ.get("BR_SWEEP", "Interpolating")
        npts = int(os.environ.get("BR_NPTS", "20001"))
        kw = dict(interpolation_tol=0.02, interpolation_max_solutions=250) \
            if kind == "Interpolating" else {}
        h.create_linear_count_sweep(setup="Setup1", unit="GHz",
                                    start_frequency=f0_ghz, stop_frequency=f1_ghz,
                                    num_of_freq_points=npts, name="Sweep1",
                                    sweep_type=kind, save_fields=False, **kw)
        R.setup_and_sweep.last_sweep = "Sweep1"
        print(f"   adapt @ {f_adapt} GHz; sweep {f0_ghz:.5f}-{f1_ghz:.5f} GHz "
              f"({npts} out pts, lw~{lw*1e3:.3f} MHz)", flush=True)
        t0 = time.time()
        h.analyze_setup("Setup1", cores=CORES)
        solve_s = time.time() - t0

        import skrf
        tpath = os.path.join(outdir, f"{tag}.s1p")
        h.export_touchstone(setup="Setup1", sweep="Sweep1", output_file=tpath)
        ntwk = skrf.Network(tpath)
        f = ntwk.f.astype(float)
        s11 = ntwk.s[:, 0, 0]

        np.savetxt(os.path.join(outdir, f"{tag}_s11.csv"),
                   np.column_stack([f, 20 * np.log10(np.abs(s11)),
                                    np.degrees(np.angle(s11)),
                                    s11.real, s11.imag]),
                   delimiter=",",
                   header="freq_hz,S11_mag_db,S11_phase_deg,S11_re,S11_im",
                   comments="")

        gd = _group_delay_qc(f, s11)
        fit = _oneport_fit(f, s11)
        res = {
            "R_ohm": Rohm, "tag": tag,
            "analytic_Qc": ANALYTIC_QC.get(Rohm),
            "analytic_formula": "Qc = (pi/4) * (R / 49.961)",
            "palace_eigenmode_Qc": PALACE_QC.get(Rohm),
            "group_delay_method": gd,
            "oneport_fit": fit,
            "si_tandelta": float(os.environ.get("SI_TANDELTA", "0")),
            "fr_ghz": round(fit.get("fr_hz", gd["fr_hz"]) / 1e9, 5),
            "Qc_group_delay": round(gd["Qc"], 1),
            "Ql_oneport_fit": round(fit.get("Ql", 0), 1),
            "Qc_oneport_fit": round(fit.get("Qc", 0), 1),
            "Qi_oneport_fit": (round(fit["Qi"], 1) if fit.get("Qi") else None),
            "n_freq": int(len(f)), "sweep_ghz": [round(f0_ghz, 5), round(f1_ghz, 5)],
            "solve_s": round(solve_s, 1),
        }
        try:
            cpath = os.path.join(outdir, f"{tag}_convergence.csv")
            h.export_convergence("Setup1", output_file=cpath)
            with open(cpath) as fp:
                res["convergence_tail"] = fp.read().strip().splitlines()[-6:]
        except Exception as e:  # noqa: BLE001
            res["convergence_error"] = repr(e)
        with open(os.path.join(outdir, f"{tag}_result.json"), "w") as fp:
            json.dump(res, fp, indent=2)
        print(json.dumps({k: v for k, v in res.items()
                          if k not in ("group_delay_method", "oneport_fit",
                                       "convergence_tail")}, indent=2), flush=True)
        return res
    finally:
        R.close_hfss(h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("R", nargs="*", help="R values (default: all three)")
    ap.add_argument("--graphical", action="store_true", default=True)
    ap.add_argument("--f0", type=float, default=7.4)
    ap.add_argument("--f1", type=float, default=8.3)
    ap.add_argument("--npts", type=int, default=8001)
    ap.add_argument("--fs", type=float, default=float(os.environ.get("BR_FS", "7.95")))
    ap.add_argument("--kind", default="Discrete")
    a = ap.parse_args()
    rs = [float(x) for x in a.R] if a.R else R_VALUES
    allres = {}
    for Rohm in rs:
        try:
            allres[f"{Rohm:.0e}"] = run_one(Rohm, HERE, a.graphical,
                                            a.f0, a.f1, a.npts, a.fs, a.kind)
        except Exception as e:  # noqa: BLE001
            allres[f"{Rohm:.0e}"] = {"error": repr(e),
                                     "traceback": traceback.format_exc()}
            print(allres[f"{Rohm:.0e}"]["traceback"], flush=True)
        with open(os.path.join(HERE, "_all_results.json"), "w") as fp:
            json.dump(allres, fp, indent=2, default=str)


if __name__ == "__main__":
    main()
