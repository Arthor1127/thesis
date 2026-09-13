"""
Per-gap linear decomposition  C(a,g) = s(g)*a + b(g).
Plot s(g) and b(g) vs g, with candidate fits and 1-sigma bands on the
intercept (from the per-gap polyfit covariance -- b is an extrapolation from
a = 60..220 um, so it carries real uncertainty).

Writes island_cap_bg.png
"""
import os, json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
rows = [json.load(open(f)) for f in sorted(glob.glob(f"{HERE}/coeff/results/*.json"))]
a = np.array([r["radius_um"] for r in rows])
g = np.array([r["gap_um"] for r in rows])
C = np.array([r["C_island_ground_fF"] for r in rows])

gv = np.array(sorted(set(g)))
S = np.zeros_like(gv, float); B = np.zeros_like(gv, float)
Serr = np.zeros_like(gv, float); Berr = np.zeros_like(gv, float)
for i, gi in enumerate(gv):
    m = g == gi
    (s, b), cov = np.polyfit(a[m], C[m], 1, cov=True)
    S[i], B[i] = s, b
    Serr[i], Berr[i] = np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1])

gg = np.linspace(gv.min(), gv.max(), 200)

# candidate fits to b(g)
cb_lin = np.polyfit(gv, B, 1)                                   # b = m g + c
A_sqrt = np.column_stack([np.ones_like(gv), np.sqrt(gv)])
cb_sqrt, *_ = np.linalg.lstsq(A_sqrt, B, rcond=None)            # b = c0 + c1 sqrt g
A_ln = np.column_stack([np.ones_like(gv), np.log(gv)])
cb_ln, *_ = np.linalg.lstsq(A_ln, B, rcond=None)                # b = c0 + c1 ln g

fig, ax = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

# ---- s(g) ----
ax[0].errorbar(gv, S, yerr=Serr, fmt="o", ms=8, capsize=4, color="C0")
cs = np.polyfit(np.log(gv), S, 1)
ax[0].plot(gg, cs[0] * np.log(gg) + cs[1], "-", color="C0",
           label=fr"$s = {cs[1]:.3f} {cs[0]:+.3f}\,\ln g$")
ax[0].set_xlabel("gap $g$ (µm)"); ax[0].set_ylabel(r"slope $s(g)$  (fF/µm)")
ax[0].set_title(r"$s(g)$: $dC/da$ at fixed gap"); ax[0].grid(alpha=.3); ax[0].legend()

# ---- b(g) ----
ax[1].errorbar(gv, B, yerr=Berr, fmt="o", ms=8, capsize=4, color="C3", label="per-gap fit intercept", zorder=5)
ax[1].plot(gg, np.polyval(cb_lin, gg), "-", color="C0",
           label=fr"linear: $b = {cb_lin[1]:.2f} + {cb_lin[0]:.3f}\,g$")
ax[1].plot(gg, cb_sqrt[0] + cb_sqrt[1] * np.sqrt(gg), "--", color="C2",
           label=fr"$\sqrt{{g}}$: $b = {cb_sqrt[0]:.2f} + {cb_sqrt[1]:.2f}\,\sqrt{{g}}$")
ax[1].plot(gg, cb_ln[0] + cb_ln[1] * np.log(gg), ":", color="C4",
           label=fr"$\ln g$: $b = {cb_ln[0]:.2f} + {cb_ln[1]:.2f}\,\ln g$")
ax[1].axhline(0, color="k", lw=.8)
ax[1].set_xlabel("gap $g$ (µm)"); ax[1].set_ylabel(r"intercept $b(g)$  (fF)")
ax[1].set_title(r"$b(g)$: $C(a,g) = s(g)\,a + b(g)$")
ax[1].grid(alpha=.3); ax[1].legend(fontsize=9)

fig.suptitle("Circular island: per-gap linear decomposition of C(a, g)   "
             r"($\varepsilon_r$=11.45, order-2, a = 60–220 µm)", fontsize=12)
fig.savefig(f"{HERE}/island_cap_bg.png", dpi=140)
print("wrote island_cap_bg.png")
print(f"\n g(um)   s(g)+-err          b(g)+-err")
for gi, s, se, b, be in zip(gv, S, Serr, B, Berr):
    print(f" {gi:5.0f}   {s:.4f} +- {se:.4f}   {b:+7.2f} +- {be:.2f}")
