#!/usr/bin/env python
"""Write params.txt: one (cpw_width_um, cpw_gap_um, tag) row per grid point,
read by submit_capacitance.sh / submit_inductance.sh via
``sed -n "${SGE_TASK_ID}p" params.txt`` (same convention as
scripts/python/island_cap/params.txt). Both submit scripts share this one file
so a given SGE_TASK_ID always means the same (width, gap) point in both sweeps.

Default grid is centred on this repo's common ~50 Ohm CPW point
(cpw_width=20um, cpw_gap=12um, e.g. notebooks/design/cap_coupler/build_design.py)
and spans it on both sides.

    python gen_params.py                       # writes params.txt, 30 points
    python gen_params.py --widths 10,20,30 --gaps 6,12,20
"""
import argparse

DEFAULT_WIDTHS_UM = [8, 12, 16, 20, 26, 34]
DEFAULT_GAPS_UM = [5, 7, 10, 14, 20]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--widths', default=None, help='comma-separated cpw_width_um list')
    ap.add_argument('--gaps', default=None, help='comma-separated cpw_gap_um list')
    ap.add_argument('--out', default='params.txt')
    args = ap.parse_args()

    widths = [float(x) for x in args.widths.split(',')] if args.widths else DEFAULT_WIDTHS_UM
    gaps = [float(x) for x in args.gaps.split(',')] if args.gaps else DEFAULT_GAPS_UM

    with open(args.out, 'w') as f:
        for w in widths:
            for g in gaps:
                tag = f"w{w:g}_g{g:g}"
                f.write(f"{w:g} {g:g} {tag}\n")

    n = len(widths) * len(gaps)
    print(f"wrote {args.out}: {len(widths)} widths x {len(gaps)} gaps = {n} points")
    print(f"widths_um = {widths}")
    print(f"gaps_um   = {gaps}")
    print(f"submit as an SGE array job:  #$ -t 1-{n}")


if __name__ == '__main__':
    main()
