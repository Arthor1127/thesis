"""Complex notch (hanger) resonator fit for a driven-S21 trace -- Phase D of the
Qc investigation (see ../QC_STATUS.md / ../HANDOFF.md).

Why this exists: qiskit_metal's `fit_transmission` (used by run_qc_sim.py) fits a
real |S21| Lorentzian and, on this design's intentionally weak coupling, was
defeated by shallow notches sitting on a sloped background. PalaceForCQED's
route instead fits the *complex* S21 with a diameter-corrected hanger model
(Khalil et al. 2012 / Probst et al. 2015), pre-removing cable delay and
environment gain/phase so the background is flattened by construction rather
than fought. This is the transmission-notch analogue of PalaceForCQED's
`notebooks/mypy.py::S11exp_no_fano_with_amp`.

Model (Probst 2015, eq. for a notch-type / side-coupled resonator):

    S21(f) = env(f) * [ 1 - (Ql/|Qc|) * e^{i*phi}
                            / (1 + 2i*Ql*(f/fr - 1)) ]
    env(f) = a * e^{i*alpha} * e^{-2*pi*i*f*tau}

Fit parameters: a, alpha, tau, fr, Ql, absQc, phi.
Derived:
    Qc = absQc / cos(phi)              (real external Q, diameter-corrected)
    1/Qi = 1/Ql - cos(phi)/absQc       (internal Q)

Pure numpy/scipy -- lmfit / resonator_tools / scikit-rf are NOT installed in the
cluster's qcg-quantum-design env.
"""
import numpy as np
from scipy.optimize import least_squares


def _env(f, a, alpha, tau):
    return a * np.exp(1j * (alpha - 2 * np.pi * f * tau))


def notch_model(f, a, alpha, tau, fr, Ql, absQc, phi):
    resonance = 1.0 - (Ql / absQc) * np.exp(1j * phi) / (1.0 + 2j * Ql * (f / fr - 1.0))
    return _env(f, a, alpha, tau) * resonance


def _estimate_delay(f, z):
    """Cable delay from the slope of the unwrapped phase over the (near-flat)
    off-resonant baseline -- use the outer 15% of points on each side."""
    n = len(f)
    k = max(3, n // 7)
    idx = np.r_[0:k, n - k:n]
    ph = np.unwrap(np.angle(z[idx]))
    slope = np.polyfit(f[idx], ph, 1)[0]
    return -slope / (2 * np.pi)


def _initial_guess(f, z):
    tau = _estimate_delay(f, z)
    z_de = z * np.exp(2j * np.pi * f * tau)
    a = float(np.median(np.abs(z_de[np.r_[0:len(f) // 7, -len(f) // 7:]])))
    if a <= 0:
        a = float(np.median(np.abs(z_de)))
    alpha = float(np.angle(np.mean(z_de[np.r_[0:len(f) // 7, -len(f) // 7:]])))
    s21_norm = z_de / (a * np.exp(1j * alpha))
    depth = 1.0 - np.abs(s21_norm)
    i0 = int(np.argmax(depth))
    fr = float(f[i0])
    d0 = float(np.clip(depth[i0], 1e-4, 0.99))          # ~ Ql/absQc on resonance (phi~0)
    # FWHM of |1 - s21_norm| ~ Lorentzian in power-ish; use half-depth crossing
    half = d0 / 2.0
    above = np.where(depth >= half)[0]
    if len(above) >= 2:
        fwhm = float(f[above[-1]] - f[above[0]])
    else:
        fwhm = (f[-1] - f[0]) / 10.0
    Ql = float(fr / max(fwhm, (f[1] - f[0]))) if fwhm > 0 else 1e4
    absQc = float(Ql / d0)
    return dict(a=a, alpha=alpha, tau=tau, fr=fr, Ql=Ql, absQc=absQc, phi=0.0)


def fit_notch(freq_hz, s21, p0=None):
    """Fit a complex S21 notch trace.

    Parameters
    ----------
    freq_hz : (N,) real array
    s21     : (N,) complex array
    p0      : optional dict overriding any initial-guess parameter

    Returns
    -------
    dict with keys: fr_hz, Ql, Qc, Qi, absQc, phi, a, alpha, tau,
        depth_db (on-resonance |S21| in dB relative to baseline),
        rmse (complex residual RMS, relative to baseline gain),
        success, message, and `model` (fitted complex array on freq_hz).
    """
    f = np.asarray(freq_hz, float)
    z = np.asarray(s21, complex)
    order = np.argsort(f)
    f, z = f[order], z[order]

    g = _initial_guess(f, z)
    if p0:
        g.update(p0)

    fr_span = f[-1] - f[0]
    # scale params to O(1) for the solver
    p_names = ["a", "alpha", "tau", "fr", "Ql", "absQc", "phi"]
    scales = np.array([g["a"] or 1.0, 1.0, max(abs(g["tau"]), 1e-9),
                       g["fr"], max(g["Ql"], 1.0), max(g["absQc"], 1.0), 1.0])
    x0 = np.array([g[k] for k in p_names]) / scales

    lb = np.array([1e-6, -np.pi, -np.inf, f[0], 1.0, 1.0, -np.pi / 2]) / scales
    ub = np.array([np.inf, np.pi, np.inf, f[-1], 1e9, 1e12, np.pi / 2]) / scales
    # fr bound: keep inside the window; widen Ql/absQc freely
    lb[3] = (f[0]) / scales[3]
    ub[3] = (f[-1]) / scales[3]

    def resid(x):
        p = x * scales
        m = notch_model(f, *p)
        r = (m - z) / (abs(g["a"]) or 1.0)
        return np.concatenate([r.real, r.imag])

    sol = least_squares(resid, x0, bounds=(lb, ub), method="trf",
                        x_scale="jac", max_nfev=20000, ftol=1e-14, xtol=1e-14)
    p = sol.x * scales
    a, alpha, tau, fr, Ql, absQc, phi = p

    Qc = absQc / np.cos(phi)
    inv_Qi = 1.0 / Ql - np.cos(phi) / absQc
    Qi = 1.0 / inv_Qi if inv_Qi > 0 else np.inf

    model = notch_model(f, *p)
    base = np.abs(_env(f, a, alpha, tau))
    s21_norm_min = np.min(np.abs(z) / base)
    depth_db = 20 * np.log10(max(s21_norm_min, 1e-12))
    rmse = float(np.sqrt(np.mean(np.abs((model - z) / (base.mean() or 1.0)) ** 2)))

    return dict(fr_hz=float(fr), Ql=float(Ql), Qc=float(Qc), Qi=float(Qi),
                absQc=float(absQc), phi=float(phi), a=float(a), alpha=float(alpha),
                tau=float(tau), depth_db=float(depth_db), rmse=rmse,
                success=bool(sol.success), message=str(sol.message),
                cost=float(sol.cost), model=model, freq_hz=f)


if __name__ == "__main__":
    # self-test on synthetic data spanning the regime the eigenmode method
    # can't reach: Qc from 5e4 to 1e6, weak coupling, clean (sim-like) trace.
    rng = np.random.default_rng(0)
    for true_Qc in (5e4, 1e5, 3e5, 1e6):
        fr = 7.5e9
        Qi = 9.2e4                       # dielectric limit for this design
        Ql = 1.0 / (1.0 / Qi + 1.0 / true_Qc)
        f = np.linspace(fr - 8e6, fr + 8e6, 1201)
        z = notch_model(f, 1.0, 0.2, 3e-9, fr, Ql, true_Qc, 0.0)
        z = z * (1 + 0.0005 * rng.standard_normal(f.shape))    # 0.05% amplitude noise
        r = fit_notch(f, z)
        err = 100 * (r["Qc"] / true_Qc - 1)
        print(f"true Qc={true_Qc:>10.0f}  fit Qc={r['Qc']:>11.0f} ({err:+5.1f}%)  "
              f"Qi={r['Qi']:>9.0f}  depth={r['depth_db']:+.2f}dB  rmse={r['rmse']:.2e}  "
              f"{'OK' if r['success'] else 'FAIL'}")
