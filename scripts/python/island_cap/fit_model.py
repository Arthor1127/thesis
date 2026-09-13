"""
Final empirical model for the circular island -> ground capacitance
(see island_capacitance_model.md).

    C(a,g) = s(g)*a + b(g)          [fF ; a,g in um ; eps_r = 11.45]
    s(g)   = s0 + s1*ln(g)
    b(g)   = b0 + b1*sqrt(g)

Joint least-squares fit to the 16-point order-2 sweep (job 1938614).
Writes island_cap_model.{json,png}.

Usage:  python fit_model.py [results_dir]   (default ./coeff/results)
"""
import os, sys, glob, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
RESDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "coeff", "results")
EPS0 = 8.8541878128e-12
EPS_R = 11.45

rows = [json.load(open(f)) for f in sorted(glob.glob(f"{RESDIR}/*.json"))]
a = np.array([r["radius_um"] for r in rows])
g = np.array([r["gap_um"] for r in rows])
C = np.array([r["C_island_ground_fF"] for r in rows])


def model(X, s0, s1, b0, b1):
    a, g = X
    return (s0 + s1 * np.log(g)) * a + (b0 + b1 * np.sqrt(g))


p, cov = curve_fit(model, (a, g), C, p0=[1.46, -0.22, -18, 1.65])
perr = np.sqrt(np.diag(cov))
pred = model((a, g), *p)
rel = (pred - C) / C
s0, s1, b0, b1 = p

report = {
    "form": "C_fF = (s0 + s1*ln(g_um))*a_um + (b0 + b1*sqrt(g_um))",
    "eps_r": EPS_R,
    "fitted_range": {"a_um": [float(a.min()), float(a.max())],
                     "g_um": [float(g.min()), float(g.max())],
                     "a_over_g": [float((a / g).min()), float((a / g).max())]},
    "coeffs": {"s0": float(s0), "s1": float(s1), "b0": float(b0), "b1": float(b1)},
    "coeffs_1sigma": {"s0": float(perr[0]), "s1": float(perr[1]),
                      "b0": float(perr[2]), "b1": float(perr[3])},
    "max_rel_err": float(np.max(np.abs(rel))),
    "rms_rel_err": float(np.sqrt(np.mean(rel ** 2))),
    "n_points": len(rows),
}
json.dump(report, open(os.path.join(HERE, "island_cap_model.json"), "w"), indent=2)
print(json.dumps(report, indent=2))

# ---- validation plot ----
fig, ax = plt.subplots(1, 3, figsize=(17, 5), constrained_layout=True)
gv = np.array(sorted(set(g)))
cols = plt.cm.viridis(np.linspace(0, .85, len(gv)))

# panel 0: C vs a, data + model, per gap
for gi, c in zip(gv, cols):
    m = g == gi
    o = np.argsort(a[m])
    ax[0].plot(a[m][o], C[m][o], "o", color=c, ms=7)
    xs = np.linspace(a.min(), a.max(), 50)
    ax[0].plot(xs, model((xs, np.full_like(xs, gi)), *p), "-", color=c, label=f"g = {gi:.0f} µm")
ax[0].set_xlabel("radius $a$ (µm)"); ax[0].set_ylabel("C (fF)")
ax[0].set_title("C(a, g): points = Palace, lines = model")
ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)

# panel 1: predicted vs simulated
ax[1].plot(C, pred, "o", ms=7)
lo, hi = C.min() * .95, C.max() * 1.05
ax[1].plot([lo, hi], [lo, hi], "k-", lw=.8)
ax[1].set_xlabel("C simulated (fF)"); ax[1].set_ylabel("C model (fF)")
ax[1].set_title(f"parity  (max {np.max(np.abs(rel))*100:.1f}%, rms {np.sqrt(np.mean(rel**2))*100:.1f}%)")
ax[1].grid(alpha=.3)

# panel 2: residuals vs a/g
sc = ax[2].scatter(a / g, rel * 100, c=g, cmap="viridis", s=80, edgecolor="k")
ax[2].axhline(0, color="k", lw=.8)
fig.colorbar(sc, ax=ax[2], label="gap g (µm)")
ax[2].set_xlabel("a / g"); ax[2].set_ylabel("model residual (%)")
ax[2].set_xscale("log"); ax[2].set_title("residuals")
ax[2].grid(alpha=.3, which="both")

fig.suptitle(r"$C(a,g) = (%.3f %+.4f\,\ln g)\,a + (%.2f %+.3f\,\sqrt{g})$   "
             r"[fF, µm, $\varepsilon_r$=11.45]" % (s0, s1, b0, b1), fontsize=12)
fig.savefig(os.path.join(HERE, "island_cap_model.png"), dpi=140)
print("\nwrote island_cap_model.{json,png}")
