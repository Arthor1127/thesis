"""Extract the eigenmode electric-field spatial distributions from Palace's
ParaView output into a compact figures/field_data_n{N}.npz per run, so the paper
field figures can be built from arrays instead of the multi-GB .vtu tree.

For each saved eigenmode it samples the FEM solution onto a fixed 2D grid in a
horizontal plane a few um above the metal, over the array region, and also a
dense 1D cut along the chain axis.

npz contents  (field_data_n{N}.npz)
-----------------------------------
    x_um, y_um            grid axes (1D),  shape (nx,), (ny,)
    z_slice_um            height of the sampling plane above the metal
    f_GHz                 (M,)  eigenfrequency of each mode
    Emag                  (M, ny, nx)  |E| = sqrt(|E_re|^2 + |E_im|^2)  [V/m]
    Ez_re, Ex_re          (M, ny, nx)  real part of the vertical / axial E
    mask                  (ny, nx) bool - grid points inside the FEM domain
    prof_x_um             (nx,)         chain-axis cut, x
    prof_Ez_re            (M, nx)       real E_z along y = 0 at z_slice
    array_x_um            (2,)  metal extent of the chain along x
    pitch_um             island pitch  (pad_length + pad_gap)

Usage:
  python export_field_data.py --run-dir eig/n40/eig_njj_40_Lj_1.2nH [--z-um 2]
"""

import argparse
import json
import os
import xml.etree.ElementTree as ET

import numpy as np
import pyvista as pv


def eig_freqs_ghz(run_dir):
    p = os.path.join(run_dir, "outputFiles", "eig.csv")
    rows = [r.strip().split(",") for r in open(p) if r.strip()]
    hdr = [h.strip().lower() for h in rows[0]]
    fcol = next(i for i, h in enumerate(hdr)
                if h.startswith("re{f") or h.startswith("f ") or h == "f")
    return [float(r[fcol]) for r in rows[1:]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--z-um", type=float, default=2.0)
    ap.add_argument("--half-len-um", type=float, default=None)
    ap.add_argument("--half-wid-um", type=float, default=80.0)
    ap.add_argument("--nx", type=int, default=600)
    ap.add_argument("--ny", type=int, default=160)
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    run_dir = os.path.abspath(args.run_dir)
    n = int(json.load(open(os.path.join(
        os.path.dirname(run_dir), "modes.json")))["n_junctions"])
    outdir = args.outdir or os.path.dirname(os.path.dirname(run_dir))  # eig/

    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from build_design import PAD_LENGTH_UM, PAD_GAP_UM, PAD_WIDTH_UM
    pitch_um = PAD_LENGTH_UM + PAD_GAP_UM
    array_half_um = 0.5 * n * pitch_um
    # island centres (local frame, chain centred on x=0): island j, j = 0..N
    island_x_um = (np.arange(n + 1) - 0.5 * n) * pitch_um

    # NOTE: Palace / SQDMetal-GMSH ParaView output is in MILLIMETRES. All the
    # geometry below is done in mm and converted to um only on the way out.
    z0 = args.z_um * 1e-3
    hx = (args.half_len_um or (array_half_um + 40.0)) * 1e-3   # array + 40um margin
    hy = args.half_wid_um * 1e-3

    pvd = os.path.join(run_dir, "outputFiles", "paraview", "eigenmode", "eigenmode.pvd")
    entries = sorted((int(float(ds.get("timestep"))), ds.get("file"))
                     for ds in ET.parse(pvd).getroot().iter("DataSet"))
    entries = [(i, f) for i, f in entries if i != 999]
    freqs = eig_freqs_ghz(run_dir)

    first = pv.read(os.path.join(run_dir, "outputFiles", "paraview",
                                 "eigenmode", entries[0][1]))
    hx = min(hx, 0.98 * first.bounds[1])   # never exceed the domain
    x = np.linspace(-hx, hx, args.nx)
    y = np.linspace(-hy, hy, args.ny)
    grid = pv.ImageData(
        dimensions=(args.nx, args.ny, 1),
        spacing=(x[1] - x[0], y[1] - y[0], 1.0),
        origin=(x[0], y[0], z0))

    M = len(entries)
    Emag = np.zeros((M, args.ny, args.nx), np.float32)
    Ez_re = np.zeros_like(Emag)
    Ex_re = np.zeros_like(Emag)
    prof_Ez = np.zeros((M, args.nx), np.float32)
    mask = None

    for mi, (ts, rel) in enumerate(entries):
        mesh = first if mi == 0 else pv.read(os.path.join(
            run_dir, "outputFiles", "paraview", "eigenmode", rel))
        s = grid.sample(mesh, pass_point_data=True)
        if mask is None:
            vm = np.asarray(s["vtkValidPointMask"]).reshape(args.ny, args.nx)
            mask = vm.astype(bool)
        Er = np.asarray(s["E_real"]).reshape(args.ny, args.nx, 3)
        Ei = np.asarray(s["E_imag"]).reshape(args.ny, args.nx, 3)
        Emag[mi] = np.sqrt((Er ** 2).sum(-1) + (Ei ** 2).sum(-1))
        Ez_re[mi] = Er[..., 2]
        Ex_re[mi] = Er[..., 0]
        prof_Ez[mi] = Er[args.ny // 2, :, 2]

    # island-averaged modal amplitude: mean real E_z over each island footprint
    # (z just above the metal) ~ island charge Q_j ~ modal coordinate. This is
    # the clean standing-wave pattern; prof_Ez_re still has the per-island
    # near-field ripple on top.
    x_um_grid = x * 1e3
    y_um_grid = y * 1e3
    hw_x = 0.35 * PAD_LENGTH_UM
    hw_y = 0.35 * PAD_WIDTH_UM
    ymask = np.abs(y_um_grid) < hw_y
    island_Ez = np.zeros((M, n + 1), np.float32)
    for j, xc in enumerate(island_x_um):
        xmask = np.abs(x_um_grid - xc) < hw_x
        if xmask.any():
            island_Ez[:, j] = Ez_re[:, ymask][:, :, xmask].mean(axis=(1, 2))

    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, f"field_data_n{n}.npz")
    np.savez_compressed(
        out,
        x_um=x * 1e3, y_um=y * 1e3, z_slice_um=args.z_um,
        f_GHz=np.array(freqs[:M]),
        Emag=Emag, Ez_re=Ez_re, Ex_re=Ex_re, mask=mask,
        prof_x_um=x * 1e3, prof_Ez_re=prof_Ez,
        island_x_um=island_x_um, island_Ez_re=island_Ez,
        array_x_um=np.array([-array_half_um, array_half_um]),
        pitch_um=pitch_um,
    )
    print(f"wrote {out}   ({M} modes, grid {args.nx}x{args.ny} @ z=+{args.z_um}um)")


if __name__ == "__main__":
    main()
