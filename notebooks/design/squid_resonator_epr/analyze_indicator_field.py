"""
Load Palace's per-element AMR error-indicator field (Cycle000003 of the
"eigenmode" ParaView output, which carries `Indicator`+`Rank` instead of
solution fields) and report both element count AND total/mean indicator
value per geometric region - the actual quantity driving refinement
decisions, not just element density.
"""
import glob
import sys

import numpy as np
import pyvista as pv

CYCLE_DIR = sys.argv[1]

PAD_RADIUS_MM = 0.090
PAD_MARGIN_MM = 0.020
SQUID_ARM_BOX_MM = (-0.710, -0.025, -0.080, 0.025)
READOUT_BOX_MM = (-1.2, 0.15, 1.2, 2.0)


def main():
    files = sorted(glob.glob(f"{CYCLE_DIR}/proc*.vtu"))
    print(f"Loading {len(files)} rank files from {CYCLE_DIR} ...")
    blocks = [pv.read(f) for f in files]
    merged = blocks[0]
    for b in blocks[1:]:
        merged = merged.merge(b)

    # Indicator/Rank are POINT data here (per-node, not per-cell)
    indicator = np.asarray(merged.point_data["Indicator"])
    centers = merged.points
    cx, cy = centers[:, 0], centers[:, 1]
    print(f"Total points: {len(indicator)}")

    x1, y1, x2, y2 = SQUID_ARM_BOX_MM
    in_squid = (cx >= x1) & (cx <= x2) & (cy >= y1) & (cy <= y2)
    in_pad = (~in_squid) & (cx * cx + cy * cy <= (PAD_RADIUS_MM + PAD_MARGIN_MM) ** 2)
    rx1, ry1, rx2, ry2 = READOUT_BOX_MM
    in_readout = (~in_squid) & (~in_pad) & (cx >= rx1) & (cx <= rx2) & (cy >= ry1) & (cy <= ry2)
    in_other = ~(in_squid | in_pad | in_readout)

    regions = {"pad": in_pad, "squid_arm": in_squid, "readout": in_readout, "other": in_other}
    total_indicator = indicator.sum()
    print(f"\n{'region':>10} {'n_points':>10} {'sum(Indicator)':>16} {'%of total ind.':>15} {'mean(Indicator)':>16}")
    for name, mask in regions.items():
        n = int(mask.sum())
        s = float(indicator[mask].sum())
        mean = s / n if n else 0.0
        pct = 100 * s / total_indicator if total_indicator else 0.0
        print(f"{name:>10} {n:>10} {s:>16.4e} {pct:>14.2f}% {mean:>16.4e}")

    # Top 1% highest-indicator points - where would AMR refine FIRST
    n_top = max(1, len(indicator) // 100)
    top_idx = np.argsort(indicator)[-n_top:]
    top_regions = {name: int(mask[top_idx].sum()) for name, mask in regions.items()}
    print(f"\nTop 1% highest-indicator points ({n_top} cells) by region: {top_regions}")


if __name__ == "__main__":
    main()
