"""Dump every data point behind figures/dispersion_fit.png into a single
figures/dispersion_data.npz (+ a flat CSV), so the paper figures can be rebuilt
without re-running anything.

npz contents
------------
scalars / fit results:
    L_J_H, C_g_F, C_Jeff_F, delta, rms_rel,
    Z_array_ohm, f_plasma_GHz, f_bandedge_GHz
per array length N in {10,23,40,100}  (prefix "n{N}_"):
    m                (int)    mode index, 1..M          (array modes only)
    f_GHz            (float)  eigenfrequency
    k_over_pi        (float)  (m + delta) / (N + 1)      open / free-free BC
    p_abs_sum        (float)  sum_j |p_j|  (~1 => array mode)
    p_per_junction   (float)  (M, N) participation matrix from port-EPR
    resid_pct        (float)  100*(f_sim - f_fit)/f_fit
    box_f_GHz        (float)  frequencies of the dropped box/package modes
dispersion fit curve (dense, model evaluated on a k grid):
    fit_k_over_pi, fit_f_GHz
"""

import csv
import json
import os

import numpy as np

from fit_dispersion import (HERE, L_J, NS, load_modes, omega_model, k_open,
                            global_fit)

OUT_NPZ = f"{HERE}/figures/dispersion_data.npz"
OUT_CSV = f"{HERE}/figures/dispersion_data.csv"


def main():
    data, box = {}, {}
    for n in NS:
        data[n], box[n] = load_modes(n)

    fit = global_fit(data, fit_delta=True)
    C_g, C_Jeff, delta, rms = fit["C_g"], fit["C_Jeff"], fit["delta"], fit["rms"]

    out = dict(
        L_J_H=L_J, C_g_F=C_g, C_Jeff_F=C_Jeff, delta=delta, rms_rel=rms,
        Z_array_ohm=np.sqrt(L_J / C_g),
        f_plasma_GHz=1.0 / (2 * np.pi * np.sqrt(L_J * C_Jeff)) / 1e9,
        f_bandedge_GHz=omega_model(np.pi, C_g, C_Jeff) / (2 * np.pi) / 1e9,
        NS=np.array(NS),
    )
    rows = []
    for n in NS:
        fs = data[n]
        M = len(fs)
        m = np.arange(1, M + 1)
        k_pi = k_open(m, n, delta) / np.pi
        f_fit = omega_model(k_open(m, n, delta), C_g, C_Jeff) / (2 * np.pi)
        resid = 100 * (fs - f_fit) / f_fit

        d = json.load(open(f"{HERE}/eig/n{n}/modes.json"))
        modes = sorted((mm for mm in d["modes"]
                        if mm.get("p_array_abs_sum", 0) > 0.5),
                       key=lambda mm: mm["f_ghz"])
        p_abs = np.array([mm["p_array_abs_sum"] for mm in modes])
        p_mat = np.array([mm.get("p_per_junction", [np.nan] * n) for mm in modes])

        out[f"n{n}_m"] = m
        out[f"n{n}_f_GHz"] = fs / 1e9
        out[f"n{n}_k_over_pi"] = k_pi
        out[f"n{n}_p_abs_sum"] = p_abs
        out[f"n{n}_p_per_junction"] = p_mat
        out[f"n{n}_resid_pct"] = resid
        out[f"n{n}_box_f_GHz"] = box[n] / 1e9

        for i in range(M):
            rows.append(dict(N=n, m=int(m[i]), f_GHz=fs[i] / 1e9,
                             k_over_pi=k_pi[i], p_abs_sum=p_abs[i],
                             resid_pct=resid[i]))

    kk = np.linspace(1e-4, np.pi, 600)
    out["fit_k_over_pi"] = kk / np.pi
    out["fit_f_GHz"] = omega_model(kk, C_g, C_Jeff) / (2 * np.pi) / 1e9

    os.makedirs(f"{HERE}/figures", exist_ok=True)
    np.savez(OUT_NPZ, **out)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["N", "m", "f_GHz", "k_over_pi",
                                          "p_abs_sum", "resid_pct"])
        w.writeheader()
        w.writerows(rows)

    print(f"wrote {OUT_NPZ}")
    print(f"wrote {OUT_CSV}  ({len(rows)} rows)")
    print("\nkeys:")
    for k in sorted(out):
        v = out[k]
        shp = getattr(v, "shape", None)
        print(f"  {k:24s} {str(shp) if shp else v}")


if __name__ == "__main__":
    main()
