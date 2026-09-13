"""Offline S11 fit for a saved .s1p -- used when the HFSS driver crashes in
post-processing but the touchstone export succeeded.  Self-contained: the two
estimators are copies of `_group_delay_qc` / `_oneport_fit` in
`run_bareres_s11.py` (which is not imported here because its module pulls in
pyaedt)."""
import json
import os
import sys

import numpy as np
import skrf
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYTIC_QC = {2e5: 3144, 1e6: 15720, 5e6: 78601}
PALACE_QC = {2e5: 3841.1, 1e6: 18495.4, 5e6: 78181.5}


def _group_delay_qc(f, s11):
    ph = np.unwrap(np.angle(s11))
    tau_g = -np.gradient(ph, f) / (2 * np.pi)
    j = int(np.argmax(np.abs(tau_g - np.median(tau_g))))
    fr = f[j]
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
    f0 = f.mean()
    mag = np.abs(s11)
    ph = np.unwrap(np.angle(s11))
    tau = -np.gradient(ph, f) / (2 * np.pi)
    j = int(np.argmax(np.abs(tau - np.median(tau))))
    fr0 = float(f[j])
    pk = abs(tau[j] - np.median(tau))
    inw = f[np.abs(tau - np.median(tau)) > 0.5 * pk]
    bw = (inw.max() - inw.min()) if inw.size > 2 else np.ptp(f) / 50
    ql0 = float(np.clip(fr0 / max(bw, f[1] - f[0]), 30, 5e6))

    def model(pp):
        fr, ql, qc_abs, phi, a, alpha, t = pp
        env = a * np.exp(1j * alpha) * np.exp(-2j * np.pi * (f - f0) * t)
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
    fr, ql, qc_abs, phi, a, alpha, t = best.x
    ql, qc_abs = abs(ql), abs(qc_abs)
    qc_complex = qc_abs * np.exp(1j * phi)
    re_inv_qc = np.real(1.0 / qc_complex)
    qi = 1.0 / (1.0 / ql - re_inv_qc) if (1.0 / ql - re_inv_qc) > 0 else None
    return dict(fr_hz=float(fr), Ql=float(ql), Qc=float(qc_abs),
                phi_rad=float(phi), Qi=(float(qi) if qi else None),
                fit_cost=float(best.cost))


tag = sys.argv[1] if len(sys.argv) > 1 else "bareres_R5e06"
Rohm = float(sys.argv[2]) if len(sys.argv) > 2 else 5e6

ntwk = skrf.Network(os.path.join(HERE, f"{tag}.s1p"))
f = ntwk.f.astype(float)
s11 = ntwk.s[:, 0, 0]

gd = _group_delay_qc(f, s11)
fit = _oneport_fit(f, s11)
res = {
    "R_ohm": Rohm, "tag": tag,
    "analytic_Qc": ANALYTIC_QC.get(Rohm),
    "analytic_formula": "Qc = (pi/4) * (R / 49.961)",
    "palace_eigenmode_Qc": PALACE_QC.get(Rohm),
    "group_delay_method": gd,
    "oneport_fit": fit,
    "si_tandelta": 1.2e-5,
    "fr_ghz": round(fit.get("fr_hz", gd["fr_hz"]) / 1e9, 5),
    "Qc_group_delay": round(gd["Qc"], 1),
    "Ql_oneport_fit": round(fit.get("Ql", 0), 1),
    "Qc_oneport_fit": round(fit.get("Qc", 0), 1),
    "Qi_oneport_fit": (round(fit["Qi"], 1) if fit.get("Qi") else None),
    "n_freq": int(len(f)),
    "sweep_ghz": [round(f[0] / 1e9, 5), round(f[-1] / 1e9, 5)],
    "solve_note": "HFSS driven interpolating sweep, fine mesh (MIN/MAXPASS "
                  "16/20, MaxDeltaS 1e-4), window centred on f_r; fit run "
                  "offline from the saved .s1p (driver crashed post-export).",
}
cpath = os.path.join(HERE, f"{tag}_convergence.csv")
if os.path.exists(cpath):
    with open(cpath) as fp:
        res["convergence_tail"] = fp.read().strip().splitlines()[-6:]

with open(os.path.join(HERE, f"{tag}_result.json"), "w") as fp:
    json.dump(res, fp, indent=2)
print(json.dumps({k: v for k, v in res.items()
                  if k not in ("group_delay_method", "convergence_tail")},
                 indent=2))
