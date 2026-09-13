#!/usr/bin/env python
"""Write params_lengthcheck.txt: (cpw_width_um, cpw_gap_um, length_um, tag) rows
for the length-linearity check -- verifying that C_total(length) and
L_total(length) really are linear in length (with a small, length-independent
intercept from the open/short end fringe), which is the premise behind
treating C_total/length and L_total/length as "the" per-unit-length C', L' at
a single length. See check_length_linearity.py for the analysis.

Two representative (width, gap) points by default: this repo's usual ~50 Ohm
point (20/12um) and the grid corner with the largest gap-to-feature ratio
(8/20um), where end-fringe contamination of a single-length estimate should be
largest -- if linearity holds even there, it should hold everywhere inside the
main sweep grid (gen_params.py).

    python gen_lengthcheck_params.py
"""
import argparse

DEFAULT_POINTS_UM = [(20.0, 12.0), (8.0, 20.0)]
DEFAULT_LENGTHS_UM = [1000.0, 3000.0, 6000.0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--points', default=None,
                    help='comma-separated w:g pairs, e.g. "20:12,8:20"')
    ap.add_argument('--lengths', default=None, help='comma-separated length_um list')
    ap.add_argument('--out', default='params_lengthcheck.txt')
    args = ap.parse_args()

    if args.points:
        points = [tuple(float(x) for x in p.split(':')) for p in args.points.split(',')]
    else:
        points = DEFAULT_POINTS_UM
    lengths = ([float(x) for x in args.lengths.split(',')] if args.lengths
              else DEFAULT_LENGTHS_UM)

    with open(args.out, 'w') as f:
        for w, g in points:
            for length in lengths:
                tag = f"w{w:g}_g{g:g}_L{length:g}"
                f.write(f"{w:g} {g:g} {length:g} {tag}\n")

    n = len(points) * len(lengths)
    print(f"wrote {args.out}: {len(points)} (w,g) points x {len(lengths)} lengths = {n} rows")
    print(f"points_um  = {points}")
    print(f"lengths_um = {lengths}")
    print(f"submit as an SGE array job:  #$ -t 1-{n}")


if __name__ == '__main__':
    main()
