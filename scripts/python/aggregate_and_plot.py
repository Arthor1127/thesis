"""
Aggregates all Stage B (driven S21) output into two plots:

  1. The ORIGINAL S21 of the resonator: the ungrounded baseline
     (ground_pos=None), S21 (dB) vs frequency, across all 4 harmonics
     stitched together.

  2. The EVOLUTION colormap: x = grounding position (mm), y = frequency
     (GHz), z = |S21| (dB), showing how the S21 spectrum changes as the
     single grounding strap is moved along the resonator. Built from the
     numeric (non-baseline) sweep positions only.

Run this after Stage B has completed for all positions/harmonics.

Usage:
    python aggregate_and_plot.py --sweep-dir ~/sweep_v2 --outdir ~/sweep_v2/plots
"""
import argparse
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt


def s21_db(s21):
    return 20 * np.log10(np.abs(s21) + 1e-30)


def load_position_dir(driven_dir):
    """Load the single s_sweep.npz file in one driven_<pos> directory.
    Returns (freqs_hz, S21) or (None, None) if not found."""
    f = os.path.join(driven_dir, "s_sweep.npz")
    if not os.path.exists(f):
        return None, None
    d = np.load(f)
    freqs = d["freq_hz"]
    s21 = d["S21"]
    order = np.argsort(freqs)
    return freqs[order], s21[order]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sweep-dir", default=os.path.expanduser("~/sweep_v2"))
    p.add_argument("--outdir", default=None)
    args = p.parse_args()
    outdir = args.outdir or os.path.join(args.sweep_dir, "plots")
    os.makedirs(outdir, exist_ok=True)

    driven_dirs = sorted(glob.glob(os.path.join(args.sweep_dir, "driven_*")))
    print(f"Found {len(driven_dirs)} driven_* directories")

    baseline_freqs, baseline_s21 = None, None
    position_data = []  # list of (pos_mm, freqs_hz, s21)

    for d in driven_dirs:
        label = os.path.basename(d).replace("driven_", "")
        freqs, s21 = load_position_dir(d)
        if freqs is None:
            print(f"  SKIP {label}: no s_sweep.npz file found")
            continue
        if label == "none":
            baseline_freqs, baseline_s21 = freqs, s21
            print(f"  Loaded baseline: {len(freqs)} points, "
                  f"{freqs[0]/1e9:.2f}-{freqs[-1]/1e9:.2f} GHz")
        else:
            pos_mm = float(label)
            position_data.append((pos_mm, freqs, s21))
            print(f"  Loaded pos={pos_mm}mm: {len(freqs)} points")

    position_data.sort(key=lambda t: t[0])

    # --- Plot 1: baseline S21 ---
    if baseline_freqs is not None:
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(baseline_freqs / 1e9, s21_db(baseline_s21), lw=1.2)
        ax.set_xlabel("Frequency (GHz)")
        ax.set_ylabel("|S21| (dB)")
        ax.set_title("Baseline (ungrounded) resonator - S21")
        for n in range(1, 5):
            ax.axvline(n * 5.0, color='gray', ls='--', lw=0.6, alpha=0.6)
        fig.tight_layout()
        out_path = os.path.join(outdir, "baseline_S21.png")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"Saved {out_path}")
    else:
        print("WARNING: no baseline data found - skipping baseline plot. "
              "Did Stage B run for ground_pos=none?")

    # --- Plot 2: evolution colormap ---
    if position_data:
        # All positions should share the same frequency grid (same harmonic
        # windows/n_points for every geometry) - verify and use the first
        # entry's grid as the common y-axis.
        ref_freqs = position_data[0][1]
        n_freq = len(ref_freqs)
        positions_mm = np.array([p for p, _, _ in position_data])
        Z = np.full((len(position_data), n_freq), np.nan)

        for i, (pos_mm, freqs, s21) in enumerate(position_data):
            if len(freqs) != n_freq or not np.allclose(freqs, ref_freqs, rtol=1e-6):
                print(f"  WARNING: pos={pos_mm}mm has a different frequency grid "
                      f"than the reference - interpolating onto the common grid.")
                s21_db_vals = np.interp(ref_freqs, freqs, s21_db(s21))
            else:
                s21_db_vals = s21_db(s21)
            Z[i, :] = s21_db_vals

        fig, ax = plt.subplots(figsize=(10, 6))
        mesh = ax.pcolormesh(positions_mm, ref_freqs / 1e9, Z.T,
                              shading='auto', cmap='viridis')
        ax.set_xlabel("Grounding strap position along resonator (mm)")
        ax.set_ylabel("Frequency (GHz)")
        ax.set_title("S21 evolution vs. grounding position")
        for n in range(1, 5):
            ax.axhline(n * 5.0, color='white', ls='--', lw=0.5, alpha=0.5)
        cbar = fig.colorbar(mesh, ax=ax)
        cbar.set_label("|S21| (dB)")
        fig.tight_layout()
        out_path = os.path.join(outdir, "S21_evolution_colormap.png")
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        print(f"Saved {out_path}")

        # Also save the raw grid, in case you want to replot/analyze later
        np.savez(os.path.join(outdir, "S21_evolution_data.npz"),
                  positions_mm=positions_mm, freq_hz=ref_freqs, S21_dB=Z)
    else:
        print("WARNING: no numeric position data found - skipping colormap.")

    print(f"\nDone. Plots in {outdir}")


if __name__ == "__main__":
    main()
