#!/usr/bin/env python
"""Plots the two analytic-model fit attempts against the trustworthy
(ceiling_flag=='ok') eigenmode Qc points, per coupling_space -- see the
QC_STATUS.md discussion following "is there an analytical model we could use
to fit the most coupled cases and extrapolate": both a naive weak-coupling
power law (Qc ~ 1/length^2) and a junction+distributed-capacitance model
(Qc ~ 1/(C0+Cm*L)^2) were tested and found not to fit (R^2=0.15 for the best
simple log-linear model, unstable/nonsensical parameters for the offset
model). This makes that visually obvious rather than just reporting numbers.

Run ON THE CLUSTER: qsub submit_plot_qc_length_fits.sh
"""
import os

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo/data root, one level up from this subfolder


def power_law_model(L, slope, intercept):
    return np.exp(intercept) * L**slope


def offset_model(L, a, b):
    return 1.0 / (a * (L + b))**2


def main():
    df = pd.read_csv(os.path.join(HERE, "qc_eigenmode.csv"))
    ok = df[df["ceiling_flag"] == "ok"].copy()

    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    import thesis_style as ts
    ts.use()

    groups = [(cs, g.sort_values("coupling_length_um")) for cs, g in ok.groupby("coupling_space_um")
             if len(g) >= 4]
    n = len(groups)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    YLIM = (3e3, 3e6)  # shared across panels so one unstable fit doesn't rescale its own panel

    fig, axs = plt.subplots(nrows, ncols, figsize=ts.figsize(2, height=3.4 * nrows), squeeze=False)
    axs_flat = axs.flatten()

    L_fine = np.linspace(90, 410, 200)
    fit_rows = []
    handles_labels = None

    for i, (cs, grp) in enumerate(groups):
        ax = axs_flat[i]
        L = grp["coupling_length_um"].values
        Qc = grp["Qc"].values
        ax.semilogy(L, Qc, "o", ms=7, color=ts.CYCLE[0], zorder=5, label="data (trustworthy)")

        # naive power law
        slope, intercept = np.polyfit(np.log(L), np.log(Qc), 1)
        ax.semilogy(L_fine, power_law_model(L_fine, slope, intercept), "--",
                   color=ts.CYCLE[1], lw=1.6, label=f"power-law fit")

        # theoretical slope=-2 reference, anchored to the data's own geometric mean
        ref_intercept = np.mean(np.log(Qc)) + 2 * np.mean(np.log(L))
        ax.semilogy(L_fine, power_law_model(L_fine, -2.0, ref_intercept), ":",
                   color="0.4", lw=1.4, label="theory (slope=-2)")

        # junction+distributed offset model
        resid_offset = None
        try:
            popt, _ = curve_fit(offset_model, L, Qc, p0=[1e-4, 50], maxfev=20000)
            ax.semilogy(L_fine, offset_model(L_fine, *popt), "-.",
                       color=ts.CYCLE[2], lw=1.6, label="junction+distributed fit")
            resid_offset = float(np.std(np.log(Qc) - np.log(offset_model(L, *popt))))
        except Exception:
            pass

        resid_power = float(np.std(np.log(Qc) - np.log(power_law_model(L, slope, intercept))))
        ax.set_title(fr"$s={cs:g}\,\mu$m (n={len(grp)}, slope={slope:+.2f})", fontsize=9, pad=6)
        ax.set_xlabel(r"$\ell$ ($\mu$m)")
        ax.set_ylabel(r"$Q_c$")
        ax.set_ylim(*YLIM)
        ax.set_xlim(70, 430)
        if handles_labels is None:
            handles_labels = ax.get_legend_handles_labels()

        fit_rows.append({"coupling_space_um": cs, "n": len(grp), "power_law_slope": slope,
                         "power_law_resid_logstd": resid_power,
                         "offset_model_resid_logstd": resid_offset})

    for j in range(len(groups), len(axs_flat)):
        axs_flat[j].axis("off")

    ts.finalize(axs_flat[:len(groups)])
    fig.legend(*handles_labels, loc="lower center", ncol=4, fontsize=8,
              bbox_to_anchor=(0.5, 0.0), frameon=False)
    fig.suptitle("Analytic model fits vs. trustworthy eigenmode Qc -- both fail "
                "(fitted slope far from the theoretical -2, poor residuals)",
                fontsize=9, y=0.995)
    fig.tight_layout(rect=[0, 0.06, 1, 0.94])
    IMAGES_DIR = os.path.join(HERE, "images")
    os.makedirs(IMAGES_DIR, exist_ok=True)
    ts.save(fig, "qc_length_fits", directory=IMAGES_DIR)
    print("wrote qc_length_fits.{png,pdf}")

    fits_df = pd.DataFrame(fit_rows)
    fits_df.to_csv(os.path.join(HERE, "qc_length_fits.csv"), index=False)
    print(fits_df.to_string(index=False))


if __name__ == "__main__":
    main()
