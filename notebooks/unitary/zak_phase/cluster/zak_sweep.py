#!/usr/bin/env python3
"""
zak_sweep.py -- one sweep: dressed-band Zak phase over a 2D (y_axis x x_axis) grid.

The colormap axes y_axis / x_axis are two of the four knobs
(omega_C, chi_C, omega_r, gamma), set in the grid JSON's config. The other two
are fixed per sweep. For every grid point, with and without the SW correction,
records:  zak_b0, zak_b1, zak_comb12, gap01, gap12.

Usage:  python zak_sweep.py --index I --grid grid_gammaomegar.json --outdir results_gammaomegar
Resumable: checkpoints one y-row at a time into <outdir>/sweep_II/_partial.npz.
"""
import argparse
import json
import os
import time

import numpy as np

import zak_kernel as K

FIELDS = ("zak_b0", "zak_b1", "zak_comb12", "gap01", "gap12")
KNOBS = ("omega_C", "chi_C", "omega_r", "gamma")


def run(index, grid_path, outdir):
    grid = json.load(open(grid_path))
    cfg = grid["config"]
    sw = next(s for s in grid["sweeps"] if s["index"] == index)
    y_axis, x_axis = cfg["y_axis"], cfg["x_axis"]

    yv = np.linspace(cfg["y_min"], cfg["y_max"], cfg["n_y"])
    xv = np.linspace(cfg["x_min"], cfg["x_max"], cfg["n_x"])
    n_y, n_x = len(yv), len(xv)
    gJ1, gJ2 = cfg["Omega_J"] * (1 + cfg["Chi_J"]), -cfg["Omega_J"] * (1 - cfg["Chi_J"])
    s = (cfg["s_A"], cfg["s_B"])
    d_ph, Nk = cfg["d_ph"], cfg["Nk"]
    fixed = {k: sw[k] for k in KNOBS if k not in (y_axis, x_axis)}

    def params_at(yval, xval):
        """-> (gc1, gc2, omega_r, gamma) for grid point (y_axis=yval, x_axis=xval)."""
        v = dict(fixed)
        v[y_axis], v[x_axis] = yval, xval
        return (0.25 * v["omega_C"] * (1 + v["chi_C"]),
                0.25 * v["omega_C"] * (1 - v["chi_C"]),
                v["omega_r"], v["gamma"])

    sdir = os.path.join(outdir, f"sweep_{index:02d}")
    os.makedirs(sdir, exist_ok=True)
    part = os.path.join(sdir, "_partial.npz")

    shape = (n_y, n_x)
    data = {f"{tag}_{f}": np.full(shape, np.nan)
            for tag in ("nosw", "sw") for f in FIELDS}
    row0 = 0
    if os.path.exists(part):
        p = np.load(part)
        for kk in data:
            data[kk] = p[kk].copy()
        row0 = int(p["_row_done"])
        print(f"[sweep {index}] resume from {y_axis} row {row0}", flush=True)

    t0 = time.time()
    for i in range(row0, n_y):
        for j, xval in enumerate(xv):
            g1, g2, wr, gam = params_at(yv[i], xv[j])
            for tag, with_sw in (("nosw", False), ("sw", True)):
                res = K.zak_point(g1, g2, gJ1, gJ2, wr, gam, d_ph, with_sw, Nk=Nk, s=s)
                for f, val in zip(FIELDS, res):
                    data[f"{tag}_{f}"][i, j] = val
        if (i + 1) % 5 == 0 or i == n_y - 1:
            np.savez(part, _row_done=i + 1, **data)
            dt = time.time() - t0
            eta = dt / (i + 1 - row0) * (n_y - i - 1)
            print(f"[sweep {index}] row {i + 1}/{n_y}  ({dt:.0f}s, ETA {eta:.0f}s)", flush=True)

    meta = dict(index=index, y_axis=y_axis, x_axis=x_axis,
                y_values=list(yv), x_values=list(xv),
                Omega_C=list(yv) if y_axis == "omega_C" else None,   # back-compat
                **fixed,
                **{k: sw[k] for k in sw if k.startswith("i_")},
                **{k: cfg[k] for k in ("Omega_J", "Chi_J", "d_ph", "Nk", "s_A", "s_B",
                                       "y_min", "y_max", "n_y", "x_min", "x_max", "n_x")},
                gamma_values=cfg["gamma_values"], omega_r_values=cfg["omega_r_values"],
                chi_C_values=cfg["chi_C_values"], omega_C_values=cfg["omega_C_values"],
                fields=list(FIELDS), wall_seconds=round(time.time() - t0, 1))
    meta = {k: v for k, v in meta.items() if v is not None}
    np.savez(os.path.join(sdir, "data.npz"), **meta, **data)
    json.dump(meta, open(os.path.join(sdir, "params.json"), "w"), indent=2)
    if os.path.exists(part):
        os.remove(part)
    print(f"[sweep {index}] done in {time.time() - t0:.0f}s -> {sdir}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--grid", default="grid_gammaomegar.json")
    ap.add_argument("--outdir", default="results_gammaomegar")
    a = ap.parse_args()
    run(a.index, a.grid, a.outdir)
