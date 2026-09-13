#!/usr/bin/env python3
"""
gaps_sweep.py -- one row of the stanene gap map: fixed B_x, all epsilon_site values.

Task i (0 <= i < mesh) computes gap_2d and gap_1d for B_vals[i] over every
E_vals[j], using the verbatim notebook kernel.  Checkpoints one column at a time
into <outdir>/row_III/_partial.npz so a killed task resumes for free.  Final
output: <outdir>/row_III.npz with B, E_vals, gap2d[mesh], gap1d[mesh].

Usage:  python gaps_sweep.py --index I --grid grid.json --outdir results_stanene
"""
import argparse
import json
import os
import time

import numpy as np

import gaps_kernel as K


def compute_point(B, e, g):
    model = K.stanene_pythtb(e, B=[B, 0.0, 0.0])
    gap2d = K.gap_2d(model, mallado=g["gap2d_mallado"])
    ribbon = model.make_finite(
        periodic_dirs=[0],
        num_cells=[g["ribbon_num_cells"]],
        glue_edges=[False],
    )
    gap1d = K.gap_1d(ribbon, mallado=g["gap1d_mallado"])
    return gap2d, gap1d


def run(index, grid_path, outdir):
    g = json.load(open(grid_path))
    mesh = g["mesh"]
    B_vals = np.linspace(g["B_min"], g["B_max"], mesh)
    E_vals = np.linspace(g["E_min"], g["E_max"], mesh)
    B = float(B_vals[index])

    sdir = os.path.join(outdir, f"row_{index:03d}")
    os.makedirs(sdir, exist_ok=True)
    part = os.path.join(sdir, "_partial.npz")

    gap2d = np.full(mesh, np.nan)
    gap1d = np.full(mesh, np.nan)
    j0 = 0
    if os.path.exists(part):
        p = np.load(part)
        gap2d, gap1d = p["gap2d"].copy(), p["gap1d"].copy()
        j0 = int(p["_col_done"])
        print(f"[row {index}] resume from col {j0}", flush=True)

    t0 = time.time()
    for j in range(j0, mesh):
        gap2d[j], gap1d[j] = compute_point(B, float(E_vals[j]), g)
        if (j + 1) % 5 == 0 or j == mesh - 1:
            np.savez(part, _col_done=j + 1, gap2d=gap2d, gap1d=gap1d)
            dt = time.time() - t0
            eta = dt / (j + 1 - j0) * (mesh - j - 1)
            print(f"[row {index}] col {j + 1}/{mesh}  ({dt:.0f}s, ETA {eta:.0f}s)",
                  flush=True)

    np.savez(os.path.join(sdir, f"../row_{index:03d}.npz"),
             index=index, B=B, E_vals=E_vals, gap2d=gap2d, gap1d=gap1d,
             wall_seconds=round(time.time() - t0, 1))
    if os.path.exists(part):
        os.remove(part)
    print(f"[row {index}] done in {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--grid", default="grid.json")
    ap.add_argument("--outdir", default="results_stanene")
    a = ap.parse_args()
    run(a.index, a.grid, a.outdir)
