#!/usr/bin/env python3
"""
make_plots.py -- two-colour trivial / topological maps of the dressed-band Zak phase.

Reads results_*/sweep_XX/data.npz (from zak_sweep.py, any pair of colormap axes)
and writes:
  * per sweep : sweep_XX.png  -- 2x2 panels (noSW/SW  x  band0/band1)
  * overview  : montage_b<band>_<tag>.png -- 8x8 grid of the sweeps

The Zak phase is quantised (0 or pi); each cell is coloured by which:
    TRIVIAL  gamma_Zak ~ 0    -> C_TRIV
    TOPO     gamma_Zak ~ pi   -> C_TOPO
Cells that are not cleanly quantised (|gamma|/pi in [0.25,0.75]) or where the
band1<->band2 gap is below GAP_TOL are drawn grey (unresolved).

--fake  : synthesise data (no compute) to check the style.
--swap  : swap the two colours.
"""
import argparse
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import numpy as np

C_TRIV = "#03AED2"       # trivial   (gamma_Zak = 0)
C_TOPO = "#D12052"       # topological (gamma_Zak = pi)
C_BAD = "#B9B9B9"        # not resolved
GAP_TOL = 1e-3

KNOBS = ("omega_C", "chi_C", "omega_r", "gamma")
AXLABEL = {"chi_C": r"$\chi_C$", "omega_r": r"$\omega_r$",
           "omega_C": r"$\Omega_C$", "gamma": r"$\gamma$"}


def _classify(z_over_pi, gapmask):
    """-> int array: 0 trivial, 1 topological, 2 unresolved."""
    a = np.abs(z_over_pi)
    out = np.full(a.shape, 2, dtype=int)
    out[a < 0.25] = 0
    out[a > 0.75] = 1
    if gapmask is not None:
        out[gapmask] = 2
    return out


_CMAP = ListedColormap([C_TRIV, C_TOPO, C_BAD])
_NORM = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], _CMAP.N)


def _axes_of(d):
    """(y_axis, x_axis, y_values, x_values) with back-compat for the omega_C runs."""
    x_axis = str(d["x_axis"])
    y_axis = str(d["y_axis"]) if "y_axis" in d else "omega_C"
    xvals = d["x_values"]
    yvals = d["y_values"] if "y_values" in d else d["Omega_C"]
    return y_axis, x_axis, yvals, xvals


def _disc_params(y_axis, x_axis):
    return [k for k in KNOBS if k not in (y_axis, x_axis)]


def _panel(ax, x, y, cls, xlabel, ylabel, title):
    ax.pcolormesh(x, y, cls, cmap=_CMAP, norm=_NORM, shading="nearest")
    if title:
        ax.set_title(title, fontsize=8)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=8)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=8)


def _swapc():
    global C_TRIV, C_TOPO, _CMAP
    C_TRIV, C_TOPO = C_TOPO, C_TRIV
    _CMAP = ListedColormap([C_TRIV, C_TOPO, C_BAD])


def per_sweep_figure(npz, p, path):
    d = np.load(npz)
    y_axis, x_axis, yarr, xarr = _axes_of(d)
    disc = _disc_params(y_axis, x_axis)
    fig, axs = plt.subplots(2, 2, figsize=(7.0, 6.4), sharex=True, sharey=True)
    for row, tag in enumerate(("nosw", "sw")):
        for col, band in enumerate((0, 1)):
            gm = (d[f"{tag}_gap12"] < GAP_TOL) if band == 1 else None
            cls = _classify(d[f"{tag}_zak_b{band}"] / np.pi, gm)
            _panel(axs[row, col], xarr, yarr, cls,
                   AXLABEL[x_axis] if row == 1 else "",
                   AXLABEL[y_axis] if col == 0 else "",
                   f"{'no SW' if tag == 'nosw' else 'SW'}  —  band {band}")
    ttl = ",   ".join(rf"{AXLABEL[k]} $= {p[k]:.3g}$" for k in disc)
    fig.suptitle(ttl, fontsize=10)
    fig.legend(handles=[Patch(fc=C_TRIV, label=r"trivial  ($\gamma_{\rm Zak}=0$)"),
                        Patch(fc=C_TOPO, label=r"topological  ($\gamma_{\rm Zak}=\pi$)"),
                        Patch(fc=C_BAD, label="unresolved / gapless")],
               loc="lower center", ncol=3, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.03))
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)


def montage(sweeps, band, tag, path):
    d0 = np.load(sweeps[0][0])
    y_axis, x_axis, _, _ = _axes_of(d0)
    rpar, cpar = _disc_params(y_axis, x_axis)          # rows, cols
    row_vals = sorted({p[rpar] for _, p in sweeps})
    col_vals = sorted({p[cpar] for _, p in sweeps})
    nr, nc = len(row_vals), len(col_vals)
    fig, axs = plt.subplots(nr, nc, figsize=(1.55 * nc, 1.5 * nr), sharex=True, sharey=True)
    axs = np.atleast_2d(axs)
    for npz, p in sweeps:
        d = np.load(npz)
        _, _, yarr, xarr = _axes_of(d)
        r, c = row_vals.index(p[rpar]), col_vals.index(p[cpar])
        gm = (d[f"{tag}_gap12"] < GAP_TOL) if band == 1 else None
        cls = _classify(d[f"{tag}_zak_b{band}"] / np.pi, gm)
        _panel(axs[r, c], xarr, yarr, cls, "", "", "")
        axs[r, c].set_xticks([]); axs[r, c].set_yticks([])
    for c, v in enumerate(col_vals):
        axs[0, c].set_title(rf"{AXLABEL[cpar]}$={v:.3g}$", fontsize=7)
    for r, v in enumerate(row_vals):
        axs[r, 0].set_ylabel(rf"{AXLABEL[rpar]}$={v:.3g}$", fontsize=7)
    fig.suptitle(rf"Zak phase, {'SW' if tag == 'sw' else 'no SW'} — band {band}"
                 rf"     (each tile: {AXLABEL[x_axis]} horiz, {AXLABEL[y_axis]} vert)", fontsize=9)
    fig.legend(handles=[Patch(fc=C_TRIV, label="trivial"),
                        Patch(fc=C_TOPO, label="topological"),
                        Patch(fc=C_BAD, label="unresolved")],
               loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)


def _fake_npz(p, cfg, path):
    y_axis, x_axis = cfg["y_axis"], cfg["x_axis"]
    yv = np.linspace(cfg["y_min"], cfg["y_max"], cfg["n_y"])
    xv = np.linspace(cfg["x_min"], cfg["x_max"], cfg["n_x"])
    YV, XV = np.meshgrid(yv, xv, indexing="ij")
    rng = np.random.default_rng(p["index"])
    v = {k: p[k] for k in _disc_params(y_axis, x_axis)}
    v[y_axis], v[x_axis] = YV, XV
    gap0 = np.sqrt(1 + np.clip(v["omega_C"], 0, 0.99)) - np.sqrt(1 - np.clip(v["omega_C"], 0, 0.99))
    topo1 = (v["omega_r"] < gap0) & (np.asarray(v["omega_C"]) > 0.05) & (np.sign(v["chi_C"]) > 0)
    topo1 = np.broadcast_to(topo1, YV.shape)
    out = {}
    for tag in ("nosw", "sw"):
        shift = 0.03 * (tag == "sw")
        z1 = np.where(v["omega_r"] < gap0 * (1 + shift), np.pi, 0.0) * (np.sign(v["chi_C"]) > 0)
        z1 = np.broadcast_to(z1, YV.shape).copy()
        z0 = np.broadcast_to(np.where(np.sign(v["chi_C"]) < 0, np.pi, 0.0), YV.shape).copy()
        out.update({f"{tag}_zak_b0": z0, f"{tag}_zak_b1": z1,
                    f"{tag}_zak_comb12": np.zeros_like(YV),
                    f"{tag}_gap01": 0.05 + 0.1 * rng.random(YV.shape),
                    f"{tag}_gap12": np.where(topo1, 2e-3 + 0.05 * rng.random(YV.shape),
                                             0.05 + 0.1 * rng.random(YV.shape))})
    np.savez(path, y_values=yv, x_values=xv, y_axis=y_axis, x_axis=x_axis, **out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results_gammaomegar")
    ap.add_argument("--grid", default="grid_gammaomegar.json")
    ap.add_argument("--plotdir", default="plots_gammaomegar")
    ap.add_argument("--swap", action="store_true")
    ap.add_argument("--fake", action="store_true")
    a = ap.parse_args()
    if a.swap:
        _swapc()
    os.makedirs(a.plotdir, exist_ok=True)
    grid = json.load(open(a.grid))

    sweeps = []
    if a.fake:
        fd = os.path.join(a.plotdir, "_fake"); os.makedirs(fd, exist_ok=True)
        for sw in grid["sweeps"]:
            npz = os.path.join(fd, f"sweep_{sw['index']:02d}.npz")
            _fake_npz(sw, grid["config"], npz)
            sweeps.append((npz, sw))
    else:
        for sd in sorted(glob.glob(os.path.join(a.results, "sweep_*"))):
            npz, pj = os.path.join(sd, "data.npz"), os.path.join(sd, "params.json")
            if os.path.exists(npz) and os.path.exists(pj):
                sweeps.append((npz, json.load(open(pj))))

    print(f"{len(sweeps)} sweeps -> {a.plotdir}")
    for npz, p in sweeps:
        per_sweep_figure(npz, p, os.path.join(a.plotdir, f"sweep_{p['index']:02d}.png"))
    for band in (0, 1):
        for tag in ("nosw", "sw"):
            montage(sweeps, band, tag, os.path.join(a.plotdir, f"montage_b{band}_{tag}.png"))
    print("done")


if __name__ == "__main__":
    main()
