#!/usr/bin/env python3
"""
make_plots.py -- assemble the row_*.npz shards into the 2D / 1D gap maps and
plot them, in the thesis house style (directives/figure-style.md).

Writes, in --outdir:
  gap_map2D.npy, gap_map1D.npy, B_vals.npy, E_vals.npy   (same names as the notebook)
  gap_map.npz                                            (everything, one file)
  gap_map_2D.{png,pdf}, gap_map_1D.{png,pdf}

--fake synthesises shard data (no compute) to check the layout.
"""
import argparse
import glob
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import thesis_style as ts
    ts.use()
    _HAVE_TS = True
except Exception as exc:                      # pragma: no cover
    print(f"[warn] thesis_style unavailable ({exc}); falling back to defaults")
    _HAVE_TS = False


def assemble(results, grid):
    mesh = grid["mesh"]
    B_vals = np.linspace(grid["B_min"], grid["B_max"], mesh)
    E_vals = np.linspace(grid["E_min"], grid["E_max"], mesh)
    gap2d = np.full((mesh, mesh), np.nan)
    gap1d = np.full((mesh, mesh), np.nan)
    found = 0
    for f in sorted(glob.glob(os.path.join(results, "row_*.npz"))):
        d = np.load(f)
        i = int(d["index"])
        gap2d[i] = d["gap2d"]
        gap1d[i] = d["gap1d"]
        found += 1
    print(f"{found}/{mesh} rows assembled")
    return B_vals, E_vals, gap2d, gap1d, found


def _map_figure(gap, B_vals, E_vals, title, path):
    fig, ax = plt.subplots(figsize=ts.figsize(1) if _HAVE_TS else (6, 5))
    im = ax.imshow(
        gap, origin="lower", aspect="auto", vmin=0,
        extent=[E_vals[0], E_vals[-1], B_vals[0], B_vals[-1]],
        cmap="thesis_cmap" if _HAVE_TS else "viridis",
    )
    fig.colorbar(im, ax=ax, label="Gap (eV)")
    ax.set_xlabel(r"$\epsilon_{\mathrm{site}}$ (eV)")
    ax.set_ylabel(r"$B_x$")
    ax.set_title(title)
    if _HAVE_TS:
        ts.finalize(ax)
        ts.save(fig, os.path.splitext(path)[0])
    else:
        fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path)


def _fake(grid):
    mesh = grid["mesh"]
    B = np.linspace(grid["B_min"], grid["B_max"], mesh)
    E = np.linspace(grid["E_min"], grid["E_max"], mesh)
    BB, EE = np.meshgrid(B, E, indexing="ij")
    Ec = 46.0 / 450.0
    g2 = np.abs(np.sqrt(np.clip(0.15 ** 2 / (0.4 ** 2 - Ec ** 2) * (EE ** 2 - Ec ** 2), 0, None)) - BB)
    g1 = np.clip(g2 - 0.02, 0, None)
    return B, E, g2, g1, mesh


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results_stanene")
    ap.add_argument("--grid", default="grid.json")
    ap.add_argument("--outdir", default="results_stanene")
    ap.add_argument("--fake", action="store_true")
    a = ap.parse_args()
    grid = json.load(open(a.grid))
    os.makedirs(a.outdir, exist_ok=True)

    if a.fake:
        B_vals, E_vals, gap2d, gap1d, _ = _fake(grid)
    else:
        B_vals, E_vals, gap2d, gap1d, found = assemble(a.results, grid)
        if found == 0:
            raise SystemExit("no row_*.npz found")

    np.save(os.path.join(a.outdir, "gap_map2D.npy"), gap2d)
    np.save(os.path.join(a.outdir, "gap_map1D.npy"), gap1d)
    np.save(os.path.join(a.outdir, "B_vals.npy"), B_vals)
    np.save(os.path.join(a.outdir, "E_vals.npy"), E_vals)
    np.savez(os.path.join(a.outdir, "gap_map.npz"),
             gap_map2D=gap2d, gap_map1D=gap1d, B_vals=B_vals, E_vals=E_vals)

    _map_figure(gap2d, B_vals, E_vals,
                "Stanene indirect gap, 2D-periodic", os.path.join(a.outdir, "gap_map_2D.png"))
    _map_figure(gap1d, B_vals, E_vals,
                "Stanene indirect gap, 1D ribbon", os.path.join(a.outdir, "gap_map_1D.png"))


if __name__ == "__main__":
    main()
