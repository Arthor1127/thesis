#!/usr/bin/env python
"""Write params_thickness.txt: (cpw_width_um, cpw_gap_um, thickness_nm, tag)
rows for the film-thickness-correction check -- DesignHelper's
compute_cpw_geometric_parameters has a finite-metal-thickness correction (Eqs.
B6-B7 in its docstring) that the main sweep's zero-thickness ("quasi-3D")
Palace geometry never exercises (see build_design.py's module docstring). This
grid instead runs FULL3D (a real metal thickness) at a few thicknesses,
against the same two representative (width, gap) points as
gen_lengthcheck_params.py, so results compare directly against
``build_design.analytic_reference(width, gap, thickness_nm=t)`` at each t.

    python gen_thickness_params.py
"""
import argparse

DEFAULT_POINTS_UM = [(20.0, 12.0), (8.0, 20.0)]
DEFAULT_THICKNESSES_NM = [20.0, 50.0, 100.0, 150.0, 200.0, 300.0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--points', default=None,
                    help='comma-separated w:g pairs, e.g. "20:12,8:20"')
    ap.add_argument('--thicknesses', default=None, help='comma-separated thickness_nm list')
    ap.add_argument('--out', default='params_thickness.txt')
    args = ap.parse_args()

    if args.points:
        points = [tuple(float(x) for x in p.split(':')) for p in args.points.split(',')]
    else:
        points = DEFAULT_POINTS_UM
    thicknesses = ([float(x) for x in args.thicknesses.split(',')] if args.thicknesses
                  else DEFAULT_THICKNESSES_NM)

    with open(args.out, 'w') as f:
        for w, g in points:
            for t in thicknesses:
                tag = f"w{w:g}_g{g:g}_t{t:g}"
                f.write(f"{w:g} {g:g} {t:g} {tag}\n")

    n = len(points) * len(thicknesses)
    print(f"wrote {args.out}: {len(points)} (w,g) points x {len(thicknesses)} thicknesses = {n} rows")
    print(f"points_um       = {points}")
    print(f"thicknesses_nm  = {thicknesses}")
    print(f"submit as an SGE array job:  #$ -t 1-{n}")


if __name__ == '__main__':
    main()
