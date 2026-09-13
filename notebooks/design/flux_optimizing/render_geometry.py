#!/usr/bin/env python
"""Render every swept geometry to a PNG tile + run sanity checks (NO simulation).

    python render_geometry.py --mode tiles  --outdir geometry_gallery --nprocs 8 \
        --shard 1 --nshards 100
    python render_geometry.py --mode sheets --outdir geometry_gallery

Each tile shows the full metal stack the way the rest of the repo does it
(``qmview``): the ground plane filled gold, its etched gaps (CPW gap, array
ground pocket, launch pads) left white, positive metal (SQUID array + the
elbow-side ground bridge + the hairpin trace) filled blue, and the 23
flux-integration apertures outlined red.

Output:
  geometry_gallery/tiles/clr<C>_seg<S>_gap<G>_semi<M>_aw<W>.png   one per grid pt
  geometry_gallery/suspect_parts/s<shard>.txt                     flagged, per shard
  geometry_gallery/_sheets/clr<C>_aw<W>.png                       contact sheets
  geometry_gallery/SUSPECT.txt                                    merged flags
"""
import argparse
import glob
import itertools
import os
import sys
import traceback
from multiprocessing import Pool

os.environ.setdefault('QISKIT_METAL_HEADLESS', '1')
os.environ.setdefault('MPLBACKEND', 'Agg')
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shapely
from shapely.ops import unary_union

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from shapely.geometry import box as _box

from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath

from build_design import (build_design, loop_aperture_polys, _design_bounds_mm,  # noqa
                          ARRAY_COMPONENT_NAME)


def poly_patches(geom, **kw):
    """matplotlib PathPatch objects for a shapely (Multi)Polygon, holes included.
    (vendored from qmview -- not installed in the cluster env)."""
    for g in getattr(geom, 'geoms', [geom]):
        if g.is_empty:
            continue
        verts, codes = [], []
        for ring in [g.exterior, *g.interiors]:
            xy = np.asarray(ring.coords)
            verts.extend(xy)
            codes.append(MplPath.MOVETO)
            codes.extend([MplPath.LINETO] * (len(xy) - 2))
            codes.append(MplPath.CLOSEPOLY)
        yield PathPatch(MplPath(verts, codes), **kw)


def ground_plane(design, min_area_mm2=1e-7):
    """Ground-plane polygons = chip sheet minus every etched feature, both the
    ``subtract=True`` polygons AND the ``subtract=True`` paths (CPW gaps), which
    ``qmview.ground_pieces`` ignores. Largest piece first; len > 1 means the
    plane is cut into disconnected plates."""
    s = design.chips.main.size
    sx, sy = float(str(s.size_x)[:-2]), float(str(s.size_y)[:-2])
    cx, cy = float(str(s.center_x)[:-2]), float(str(s.center_y)[:-2])
    sheet = _box(cx - sx / 2, cy - sy / 2, cx + sx / 2, cy + sy / 2)

    cuts = []
    poly = design.qgeometry.tables['poly']
    cuts += list(poly[poly['subtract'].astype(bool)]['geometry'])
    path = design.qgeometry.tables['path']
    for _, r in path[path['subtract'].astype(bool)].iterrows():
        w = float(r.get('width') or 0.0) or 0.01
        cuts.append(r['geometry'].buffer(w / 2.0, cap_style=2, join_style=2))
    g = sheet.difference(unary_union(cuts)) if cuts else sheet
    return sorted((x for x in getattr(g, 'geoms', [g]) if x.area > min_area_mm2),
                  key=lambda x: -x.area)

AXES = {
    'route_clearance_um': [30., 33., 36., 40.],
    'segment_length_um':  [0., 3., 6., 9., 12.],
    'squid_gap_um':       list(np.round(np.linspace(14., 40., 6), 1)),
    'arm_semi_length_um': [8., 12., 16., 20., 24.],
    'arm_width_um':       [3., 4.5, 6., 8.],
}
KEYS = list(AXES)

GAPS = [f'{v:g}' for v in AXES['squid_gap_um']]       # sheet columns
SEMIS = [f'{v:g}' for v in AXES['arm_semi_length_um']]  # sheet row-blocks
SEGS = [f'{v:g}' for v in AXES['segment_length_um']]    # rows within a block
CLRS = [f'{v:g}' for v in AXES['route_clearance_um']]
AWS = [f'{v:g}' for v in AXES['arm_width_um']]

GROUND_FC = '#d9c37a'
GROUND_EC = '#a98f3d'
METAL_FC = '#3b7dd8'
METAL_EC = '#16375e'
APERT_C = '#c0392b'


def _tag(p):
    return (f"clr{p['route_clearance_um']:g}_seg{p['segment_length_um']:g}"
            f"_gap{p['squid_gap_um']:g}_semi{p['arm_semi_length_um']:g}"
            f"_aw{p['arm_width_um']:g}")


def _polys(tbl, cid, subtract=None):
    m = tbl['component'] == cid
    if subtract is not None:
        m &= tbl['subtract'].astype(bool) == subtract
    return list(tbl[m].geometry)


def render_one(args):
    outdir, combo = args
    p = dict(zip(KEYS, [float(x) for x in combo]))
    tag = _tag(p)
    png = os.path.join(outdir, 'tiles', tag + '.png')
    os.makedirs(os.path.dirname(png), exist_ok=True)
    flags = []
    try:
        d = build_design(chip_margin_um=250, chip_size_z_um=200,
                         include_array_metal=True, ground_near_end=True,
                         route_clearance_um=p['route_clearance_um'],
                         segment_length_um=p['segment_length_um'],
                         squid_gap_um=p['squid_gap_um'],
                         arm_semi_length_um=p['arm_semi_length_um'],
                         arm_width_um=p['arm_width_um'])
        d.rebuild()
        poly = d.qgeometry.tables['poly']
        path = d.qgeometry.tables['path']
        sub = poly['subtract'].astype(bool)
        comps = d.components

        gpieces = ground_plane(d)
        aid = comps[ARRAY_COMPONENT_NAME].id
        arr = unary_union(_polys(poly, aid, subtract=False))
        pockets = unary_union(_polys(poly, aid, subtract=True))
        route = comps['readout_1']
        rpaths = path[path['component'] == route.id]
        rline = unary_union(list(rpaths.geometry))
        try:
            rw = float(rpaths['width'].iloc[0])
        except Exception:
            rw = 0.020
        rtrace = rline.buffer(rw / 2.0, cap_style=2, join_style=2)
        try:
            bid = comps[ARRAY_COMPONENT_NAME + '_gnd_bridge'].id
            bridge = unary_union(_polys(poly, bid))
        except (KeyError, AttributeError):
            bridge = None

        apert = loop_aperture_polys(squid_gap_um=p['squid_gap_um'],
                                    arm_semi_length_um=p['arm_semi_length_um'],
                                    arm_width_um=p['arm_width_um'],
                                    segment_length_um=p['segment_length_um'])
        apert_mm = shapely.affinity.scale(apert, 1e3, 1e3, origin=(0, 0))

        # ---- sanity checks ----
        if len(apert_mm.geoms) != 23:
            flags.append(f'apertures={len(apert_mm.geoms)} (!=23)')
        if min(g.area for g in apert_mm.geoms) < 1e-8:
            flags.append('degenerate aperture')
        rl = route.length * 1e3
        if not (400 < rl < 4500):
            flags.append(f'route length {rl:.0f}um out of [400,4500]')
        if rtrace.intersection(arr).area > 5e-5:
            flags.append('hairpin trace touches array metal')
        if bridge is not None:
            if not bridge.intersects(arr):
                flags.append('ground bridge does not touch array')
            if not any(bridge.intersects(g) for g in gpieces):
                flags.append('ground bridge does not reach the ground plane')
        else:
            flags.append('ground bridge missing')
        if len(gpieces) != 1:
            flags.append(f'ground plane in {len(gpieces)} pieces (want 1)')
        ax0, ay0, ax1, ay1 = arr.bounds
        if (ax1 - ax0) < 0.3 or (ay1 - ay0) < 0.01 or (ay1 - ay0) > 0.25:
            flags.append(f'array bbox odd: {ax1-ax0:.3f} x {ay1-ay0:.3f} mm')
        for g in (arr, rtrace):
            if not np.isfinite(np.asarray(g.bounds)).all():
                flags.append('non-finite coords')

        # ---- draw ----
        fig, ax = plt.subplots(figsize=(12, 2.9))
        for piece in gpieces:                       # ground plane (gold)
            for pp in poly_patches(piece, fc=GROUND_FC, ec=GROUND_EC, lw=0.5, zorder=0):
                ax.add_patch(pp)
        for _, row in poly[~sub].iterrows():        # positive metal (blue)
            for pp in poly_patches(row['geometry'], fc=METAL_FC, ec=METAL_EC,
                                   lw=0.5, zorder=2):
                ax.add_patch(pp)
        for pp in poly_patches(rtrace, fc=METAL_FC, ec=METAL_EC, lw=0.4, zorder=2):
            ax.add_patch(pp)                        # hairpin trace metal (blue)
        for g in apert_mm.geoms:                    # flux apertures (red)
            ax.plot(*g.exterior.xy, color=APERT_C, lw=0.8, zorder=4)

        mnx, mny, mxx, mxy = _design_bounds_mm(d)
        padx, pady = (mxx - mnx) * 0.02, (mxy - mny) * 0.12
        ax.set_xlim(mnx - padx, mxx + padx)
        ax.set_ylim(mny - pady, mxy + pady)
        ax.set_aspect('equal')
        ax.set_xlabel('x [mm]', fontsize=7)
        ax.set_ylabel('y [mm]', fontsize=7)
        ax.tick_params(labelsize=6)
        ttl = (f"clr {p['route_clearance_um']:g}   seg {p['segment_length_um']:g}   "
               f"gap {p['squid_gap_um']:g}   semi {p['arm_semi_length_um']:g}   "
               f"aw {p['arm_width_um']:g}    ·   route {rl:.0f} um   ·   "
               f"gnd pieces {len(gpieces)}")
        ax.set_title(ttl + ('    !! ' + '; '.join(flags) if flags else ''),
                     fontsize=8, color=APERT_C if flags else '#1a1a2e')
        fig.savefig(png, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)
    except Exception as e:
        flags.append(f'BUILD FAILED: {type(e).__name__}: {e}')
        traceback.print_exc()
    return (tag, flags)


def do_tiles(a):
    allc = list(itertools.product(*[AXES[k] for k in KEYS]))
    combos = allc[a.shard - 1::a.nshards]
    print(f"shard {a.shard}/{a.nshards}: {len(combos)} geometries, {a.nprocs} procs",
          flush=True)
    with Pool(a.nprocs) as pool:
        res = pool.map(render_one, [(a.outdir, c) for c in combos], chunksize=4)
    sus = [(t, f) for t, f in res if f]
    sd = os.path.join(a.outdir, 'suspect_parts')
    os.makedirs(sd, exist_ok=True)
    with open(os.path.join(sd, f's{a.shard:04d}.txt'), 'w') as fh:
        for t, f in sorted(sus):
            fh.write(f"{t}\n    {'; '.join(f)}\n")
    print(f"shard {a.shard}/{a.nshards}: {len(combos)} geoms, {len(sus)} flagged",
          flush=True)


def do_sheets(a):
    """Compose one contact sheet per (route_clearance, arm_width).
    Columns = squid_gap; rows = arm_semi_length blocks x segment_length."""
    from PIL import Image
    tiles = os.path.join(a.outdir, 'tiles')
    out = os.path.join(a.outdir, '_sheets')
    os.makedirs(out, exist_ok=True)
    # tile pixel box on the sheet (keep the source aspect ~ 12:2.9)
    TW, TH = a.tile_w, int(round(a.tile_w * 2.9 / 12)) + 26
    for clr in CLRS:
        for aw in AWS:
            cols, rows = len(GAPS), len(SEMIS) * len(SEGS)
            sheet = Image.new('RGB', (cols * TW, rows * TH + 46), 'white')
            try:
                from PIL import ImageDraw
                ImageDraw.Draw(sheet).text(
                    (10, 14), f'route_clearance = {clr} um     arm_width = {aw} um'
                    f'      (cols: squid_gap {GAPS[0]}..{GAPS[-1]} um   |   '
                    f'row blocks: arm_semi_length {SEMIS[0]}..{SEMIS[-1]} um   |   '
                    f'rows in block: segment_length {SEGS[0]}..{SEGS[-1]} um)',
                    fill='black')
            except Exception:
                pass
            found = 0
            for si, semi in enumerate(SEMIS):
                for gi, seg in enumerate(SEGS):
                    r = si * len(SEGS) + gi
                    for ci, gap in enumerate(GAPS):
                        f = os.path.join(
                            tiles, f'clr{clr}_seg{seg}_gap{gap}_semi{semi}_aw{aw}.png')
                        if not os.path.exists(f):
                            continue
                        found += 1
                        im = Image.open(f).convert('RGB')
                        im.thumbnail((TW - 4, TH - 4))
                        sheet.paste(im, (ci * TW + 2, r * TH + 2 + 46))
            sheet.save(os.path.join(out, f'clr{clr}_aw{aw}.png'))
            print(f'  clr{clr}_aw{aw}: {found}/{cols*rows} tiles -> {sheet.size}',
                  flush=True)
    # merge suspect parts
    parts = sorted(glob.glob(os.path.join(a.outdir, 'suspect_parts', 's*.txt')))
    with open(os.path.join(a.outdir, 'SUSPECT.txt'), 'w') as out_fh:
        n = 0
        for pth in parts:
            txt = open(pth).read()
            n += txt.count('\n    ')
            out_fh.write(txt)
        out_fh.write(f'\n# {n} flagged geometries across {len(parts)} shards\n')
    print(f'SUSPECT.txt: {n} flagged', flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['tiles', 'sheets'], default='tiles')
    ap.add_argument('--outdir', default='geometry_gallery')
    ap.add_argument('--nprocs', type=int, default=2)
    ap.add_argument('--shard', type=int, default=1)
    ap.add_argument('--nshards', type=int, default=1)
    ap.add_argument('--tile-w', type=int, default=1400)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    if a.mode == 'tiles':
        do_tiles(a)
    else:
        do_sheets(a)


if __name__ == '__main__':
    main()
