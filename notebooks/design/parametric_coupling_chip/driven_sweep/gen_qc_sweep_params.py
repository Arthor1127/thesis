#!/usr/bin/env python
"""Write params_qc.txt: (coupling_space_um, coupling_length_um, tag) rows for
the Qc sweep, read by submit_qc_sweep.sh via
``sed -n "${SGE_TASK_ID}p" params_qc.txt`` (same convention as
notebooks/design/cpw_impedance_test/gen_params.py).

Default grid matches readout_resonators.ipynb's own geometry-sanity sweep
(coupling_length_um up to 400um, safely under COUPLING_BUDGET_UM=550um so
down_length_um never goes negative -- see readout_resonators.py).

    python gen_qc_sweep_params.py                      # writes params_qc.txt, 35 points
    python gen_qc_sweep_params.py --spaces 10,20,30 --lengths 150,250,400
"""
import argparse
import os

DEFAULT_SPACES_UM = [10.0, 15.0, 20.0, 25.0, 30.0]
DEFAULT_LENGTHS_UM = [100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0]
# writes next to this script (driven_sweep/params_qc.txt) regardless of cwd
_DEFAULT_OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'params_qc.txt')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--spaces', default=None, help='comma-separated coupling_space_um list')
    ap.add_argument('--lengths', default=None, help='comma-separated coupling_length_um list')
    ap.add_argument('--out', default=_DEFAULT_OUT)
    args = ap.parse_args()

    spaces = [float(x) for x in args.spaces.split(',')] if args.spaces else DEFAULT_SPACES_UM
    lengths = [float(x) for x in args.lengths.split(',')] if args.lengths else DEFAULT_LENGTHS_UM

    with open(args.out, 'w') as f:
        for cs in spaces:
            for cl in lengths:
                tag = f"cs{cs:g}_cl{cl:g}"
                f.write(f"{cs:g} {cl:g} {tag}\n")

    n = len(spaces) * len(lengths)
    print(f"wrote {args.out}: {len(spaces)} spaces x {len(lengths)} lengths = {n} points")
    print(f"coupling_space_um  = {spaces}")
    print(f"coupling_length_um = {lengths}")
    print(f"submit as an SGE array job:  #$ -t 1-{n}")


if __name__ == '__main__':
    main()
