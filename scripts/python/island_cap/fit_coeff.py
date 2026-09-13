"""
Refine the circular-pad capacitance model from the ORDER-2 (a,g) sweep
(job 1938614, 16 points: a in {60,100,150,220} um x g in {6,12,24,45} um,
solver order 2, res_min = max(g/6, 1.2) um, quasi-3D, eps_r = 11.9,
substrate thickness h = 500 um).

Base model:   C / (eps_eff * a)  =  A  +  B * ln(a/g)
              eps_eff = eps0 * (1 + eps_r) / 2

The naive physical form has A = 8 (isolated disc, C = 8 eps a) and
B = 8/pi ~ 2.546 (leading coplanar-gap log enhancement).  Fitting A, B over
the full a/g = 1.3 .. 37 range leaves a *systematic* ~14% residual that grows
with the island radius at fixed gap: C/(eps_eff a) is NOT a pure function of
a/g because the 500 um substrate and the finite chip box introduce their own
length scale.  Adding a linear term in a/h (h = substrate thickness) removes
most of it and lets A, B be read off in the a/h -> 0 (thick-substrate) limit:

              C / (eps_eff * a)  =  A  +  B * ln(a/g)  -  F * (a/h)

Usage:  python fit_coeff.py [results_dir]   (default ./coeff/results)
Writes: island_cap_coeff.csv, island_cap_coeff.png, island_cap_coeff.json
"""
import os
import sys
import glob
import json
import csv

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RESDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "coeff", "results")
EPS0 = 8.8541878128e-12
H_SUB = 500e-6           # substrate thickness used in run_cap_island.py (size_z)

rows = []
for fn in sorted(glob.glob(f"{RESDIR}/*.json")):
    d = json.load(open(fn)); d.pop("maxwell_matrix_F", None)
    rows.append(d)
if not rows:
    sys.exit(f"no result json in {RESDIR}")

eps_r = rows[0]["eps_r"]
eps_eff = EPS0 * (1.0 + eps_r) / 2.0

a = np.array([r["radius_um"] for r in rows]) * 1e-6          # m
g = np.array([r["gap_um"] for r in rows]) * 1e-6             # m
C = np.array([r["C_island_ground_fF"] for r in rows]) * 1e-15  # F
y = C / (eps_eff * a)                                        # dimensionless
L = np.log(a / g)


def linfit(cols, mask=None):
    m = np.ones(len(y), bool) if mask is None else mask
    M = np.column_stack(cols)[m]
    coef, *_ = np.linalg.lstsq(M, y[m], rcond=None)
    yh = M @ coef
    r2 = 1.0 - np.sum((y[m] - yh) ** 2) / np.sum((y[m] - y[m].mean()) ** 2)
    maxr = float(np.max(np.abs(yh - y[m]) / y[m]))
    rms = float(np.sqrt(np.mean(((yh - y[m]) / y[m]) ** 2)))
    return coef, dict(R2=float(r2), max_rel_resid=maxr, rms_rel_resid=rms, n=int(m.sum()))


ones = np.ones_like(L)

# fit 1: full range, plain 2-parameter log law
(A1, B1), s1 = linfit([ones, L])

# fit 2: full range, + substrate term  y = A + B ln(a/g) - F (a/h)
(A2, B2, negF2), s2 = linfit([ones, L, a / H_SUB])
F2 = -negF2

# fit 3: transmon-relevant sub-box (a >= 100 um, g <= 24 um), plain log law
box = (a >= 100e-6) & (g <= 24e-6)
(A3, B3), s3 = linfit([ones, L], mask=box)

# naive physical form, evaluated on the full set
naive = 8.0 + (8.0 / np.pi) * L
s_naive = dict(max_rel_resid=float(np.max(np.abs(naive - y) / y)),
               rms_rel_resid=float(np.sqrt(np.mean(((naive - y) / y) ** 2))))

report = {
    "source": "job 1938614, 16 pts, order 2, eps_r=%.1f, h_substrate=%g um" % (eps_r, H_SUB * 1e6),
    "eps_eff_over_eps0": eps_eff / EPS0,
    "naive_physical": {"A": 8.0, "B": 8.0 / np.pi, **s_naive},
    "fit1_plain_log_full_range": {"A": float(A1), "B": float(B1), **s1},
    "fit2_with_substrate_term": {
        "form": "C/(eps_eff a) = A + B ln(a/g) - F (a/h),  h = %g um" % (H_SUB * 1e6),
        "A": float(A2), "B": float(B2), "F": float(F2), **s2},
    "fit3_plain_log_transmon_box_a>=100_g<=24": {"A": float(A3), "B": float(B3), **s3},
    "note": ("Plain a/g scaling has a ~14 pct systematic residual that grows "
             "with a at fixed g -- the finite 500um substrate breaks the scale "
             "invariance. The a/h term absorbs it (max resid 14 pct -> 7 pct); "
             "A2,B2 are the thick-substrate coefficients. For real transmon "
             "islands (a=100-220um, g=6-24um on 500um Si) use fit3 "
             f"(A={A3:.2f}, B={B3:.2f}, 4 pct max)."),
}

# --- CSV ---
yhat1 = A1 + B1 * L
yhat2 = A2 + B2 * L - F2 * (a / H_SUB)
with open(os.path.join(HERE, "island_cap_coeff.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["radius_um", "gap_um", "a_over_g", "C_fF", "C_over_epsEff_a",
                "fit1_pred", "fit1_resid_%", "fit2_pred", "fit2_resid_%", "naive_pred", "naive_resid_%"])
    for i in np.argsort(a / g):
        w.writerow([f"{a[i]*1e6:.0f}", f"{g[i]*1e6:.0f}", f"{a[i]/g[i]:.2f}",
                    f"{C[i]*1e15:.3f}", f"{y[i]:.4f}",
                    f"{yhat1[i]:.4f}", f"{(yhat1[i]-y[i])/y[i]*100:+.1f}",
                    f"{yhat2[i]:.4f}", f"{(yhat2[i]-y[i])/y[i]*100:+.1f}",
                    f"{naive[i]:.4f}", f"{(naive[i]-y[i])/y[i]*100:+.1f}"])

# --- plot ---
fig, ax = plt.subplots(1, 3, figsize=(18, 5.2), constrained_layout=True)
gaps = sorted(set(g))
cmap = plt.cm.viridis(np.linspace(0, 0.9, len(gaps)))

# panel 0: dimensionless C vs ln(a/g), colour = gap
for gg, cc in zip(gaps, cmap):
    m = g == gg
    o = np.argsort(L[m])
    ax[0].plot(L[m][o], y[m][o], "o-", color=cc, ms=7, label=f"g = {gg*1e6:.0f} µm")
Ls = np.linspace(L.min(), L.max(), 100)
ax[0].plot(Ls, A1 + B1 * Ls, "k-", lw=2, label=fr"fit: ${A1:.2f} + {B1:.2f}\ln(a/g)$")
ax[0].plot(Ls, 8.0 + (8.0 / np.pi) * Ls, "--", color="gray",
           label=r"naive $8 + \frac{8}{\pi}\ln(a/g)$")
ax[0].set_xlabel(r"$\ln(a/g)$")
ax[0].set_ylabel(r"$C\,/\,(\varepsilon_{\rm eff}\,a)$")
ax[0].set_title(f"plain log law: $R^2$={s1['R2']:.3f}, max resid {s1['max_rel_resid']*100:.0f}%\n"
                "(spread by gap = broken $a/g$ scaling)")
ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)

# panel 1: the substrate-scale term. residual of plain fit vs a/h
resid1 = (yhat1 - y) / y * 100
for gg, cc in zip(gaps, cmap):
    m = g == gg
    o = np.argsort(a[m])
    ax[1].plot((a[m] / H_SUB)[o], resid1[m][o], "o-", color=cc, ms=7, label=f"g = {gg*1e6:.0f} µm")
ax[1].axhline(0, color="k", lw=.8)
ax[1].set_xlabel(r"$a / h$   ($h$ = substrate thickness, 500 µm)")
ax[1].set_ylabel("plain-fit residual (%)")
ax[1].set_title("residual grows with $a/h$:\nfinite substrate breaks scale invariance")
ax[1].grid(alpha=.3); ax[1].legend(fontsize=8)

# panel 2: predicted vs simulated, all three models
ax[2].plot([y.min(), y.max()], [y.min(), y.max()], "k-", lw=.8)
ax[2].plot(y, naive, "s", ms=6, color="gray", alpha=.6,
           label=f"naive (max {s_naive['max_rel_resid']*100:.0f}%)")
ax[2].plot(y, yhat1, "o", ms=6, color="C0",
           label=f"$A+B\\ln(a/g)$ (max {s1['max_rel_resid']*100:.0f}%)")
ax[2].plot(y, yhat2, "^", ms=6, color="C3",
           label=f"$+\\,F(a/h)$ (max {s2['max_rel_resid']*100:.0f}%)")
ax[2].set_xlabel(r"simulated  $C/(\varepsilon_{\rm eff} a)$")
ax[2].set_ylabel(r"predicted")
ax[2].set_title(f"substrate term: $A$={A2:.2f}, $B$={B2:.2f}, $F$={F2:.2f}\n"
                f"transmon box $a{{\\geq}}100,g{{\\leq}}24$: $A$={A3:.2f}, $B$={B3:.2f} ({s3['max_rel_resid']*100:.0f}%)")
ax[2].grid(alpha=.3); ax[2].legend(fontsize=8)

fig.suptitle(f"Circular-pad island$\\to$ground capacitance: coefficient refinement  "
             f"($\\varepsilon_r$={eps_r}, Palace order-2, {len(rows)} points)", fontsize=13)
fig.savefig(os.path.join(HERE, "island_cap_coeff.png"), dpi=140)

json.dump(report, open(os.path.join(HERE, "island_cap_coeff.json"), "w"), indent=2)
print(json.dumps(report, indent=2))
print("\nwrote island_cap_coeff.{csv,png,json}")
