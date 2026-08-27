"""
Sweep planner for the quarter-wave strap study.

WHY THIS EXISTS: a strap at interior position x0 shorts the line there,
splitting it into two independent resonators -

    A: [0, x0]   short-short  ->  f = n * v / (2*x0)
    B: [x0, L]   short-open   ->  f = (2m-1) * v / (4*(L-x0))

so the spectrum is the union of two ladders whose spacings both depend
on x0. The LOWEST mode therefore moves a lot across the sweep, and NOT
monotonically - for L=5.92905mm it runs 5.46 GHz (strap near the short
end) up to ~14 GHz around x0/L~0.7, then back down to 10.9 GHz near the
open end.

That matters because Palace's eigensolver is shift-and-invert: it
searches around config Solver.Eigenmode.Target (SQDMetal's
'starting_freq'). A single fixed target - say 5 GHz, correct only for
the ungrounded baseline - would sit far from the true mode at most strap
positions. On the half-wave design that exact mistake produced an
eig.csv with a header and ZERO converged modes.

So each array task gets its own --starting-freq and --target-freq-ghz,
predicted analytically here.

CAVEAT: the two-segment model assumes an IDEAL hard short at x0. A real
10um strap has finite width and inductance, so treat these as search
CENTRES, not expected answers - the point is to put the shift close
enough that the solver converges, then read the actual frequency off
eig.csv. Expect real modes to sit somewhat below the ideal prediction.

Usage:
    python qw_sweep_plan.py                 # print the whole plan
    python qw_sweep_plan.py --task-id 3     # emit one task's params
"""
import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design_quarterwave import (
    TOTAL_LENGTH_MM, EDGE_MARGIN_MM, sweep_positions_mm)

# Phase velocity from design_helper's elliptic-integral CPW model for
# this design's cross-section (20um trace / 12.25um gap, eps_r=11.9,
# alpha_inductance=0, t=100nm). Same number the length calibration in
# build_design_quarterwave.py is derived from - keep them consistent.
V_PHASE_M_PER_S = 1.1858e8


def predict_modes(x0_mm, total_length_mm=TOTAL_LENGTH_MM, n_modes=4):
    """Lowest few predicted mode frequencies (Hz) for a strap at x0_mm,
    measured from the SHORT end. Returns a sorted list."""
    v = V_PHASE_M_PER_S
    L = total_length_mm * 1e-3
    x0 = x0_mm * 1e-3
    if not (0 < x0 < L):
        raise ValueError(f"x0_mm={x0_mm} outside (0, {total_length_mm})")

    freqs = []
    # Segment A: short-short, full-wave ladder n*v/(2*x0)
    for n in range(1, n_modes + 1):
        freqs.append(n * v / (2 * x0))
    # Segment B: short-open, odd-quarter-wave ladder (2m-1)*v/(4*(L-x0))
    for m in range(1, n_modes + 1):
        freqs.append((2 * m - 1) * v / (4 * (L - x0)))
    return sorted(freqs)


def baseline_modes(total_length_mm=TOTAL_LENGTH_MM, n_modes=4):
    """Ungrounded quarter-wave ladder: (2n-1)*v/(4L)."""
    v = V_PHASE_M_PER_S
    L = total_length_mm * 1e-3
    return [(2 * n - 1) * v / (4 * L) for n in range(1, n_modes + 1)]


def build_plan(n_points=9, total_length_mm=TOTAL_LENGTH_MM):
    """Task list. Task 1 is always the ungrounded baseline, matching the
    half-wave pipeline's convention."""
    plan = [{
        "task_id": 1,
        "ground_pos": "none",
        "target_hz": baseline_modes(total_length_mm)[0],
    }]
    for i, x0 in enumerate(sweep_positions_mm(
            n_points=n_points, total_length_mm=total_length_mm)):
        plan.append({
            "task_id": i + 2,
            "ground_pos": x0,
            "target_hz": predict_modes(x0, total_length_mm)[0],
        })
    return plan


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-points", type=int, default=9)
    p.add_argument("--total-length-mm", type=float, default=TOTAL_LENGTH_MM)
    p.add_argument("--task-id", type=int, default=None,
                    help="If given, print just this task's params as "
                         "'<ground_pos> <starting_freq_hz> <target_ghz>' for "
                         "the array job to consume. Otherwise print the plan.")
    args = p.parse_args()

    plan = build_plan(args.n_points, args.total_length_mm)

    if args.task_id is not None:
        if not (1 <= args.task_id <= len(plan)):
            print(f"task-id {args.task_id} out of range 1..{len(plan)}",
                  file=sys.stderr)
            sys.exit(1)
        t = plan[args.task_id - 1]
        print(f"{t['ground_pos']} {t['target_hz']:.6e} "
              f"{t['target_hz']/1e9:.6f}")
        return

    L = args.total_length_mm
    print(f"L = {L} mm,  v = {V_PHASE_M_PER_S:.4e} m/s")
    print(f"Ungrounded ladder (GHz): "
          f"{[round(f/1e9, 3) for f in baseline_modes(L)]}")
    print(f"\nSpecial point 2L/3 = {2*L/3:.4f} mm - both segments land on "
          f"3*f1 there (degenerate); it is the 3*lambda/4 mode's node, so "
          f"the strap barely perturbs that mode.\n")
    print(f"{'task':>4} {'x0(mm)':>8} {'x0/L':>6} {'target(GHz)':>12}  "
          f"predicted lowest modes (GHz)")
    print("-" * 78)
    for t in plan:
        if t["ground_pos"] == "none":
            modes = baseline_modes(L)
            xs, ratio = "baseline", ""
        else:
            modes = predict_modes(t["ground_pos"], L)
            xs = f"{t['ground_pos']:.3f}"
            ratio = f"{t['ground_pos']/L:.3f}"
        print(f"{t['task_id']:>4} {xs:>8} {ratio:>6} "
              f"{t['target_hz']/1e9:>12.3f}  "
              f"{[round(f/1e9, 2) for f in modes[:4]]}")
    print(f"\nTotal tasks: {len(plan)}  ->  use '#$ -t 1-{len(plan)}'")


if __name__ == "__main__":
    main()
