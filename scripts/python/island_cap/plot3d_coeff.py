"""
3D view of the island->ground capacitance over the (a, g) grid
(job 1938614, 16 pts, Palace order-2).

NOTE on the sim setup (from the generated palace config):
  * substrate permittivity is 11.45 (not 11.9)
  * ALL SIX faces of the bounding box are grounded PEC -- including the box
    top and bottom, which sit +-500 um from the island plane. So a fixed
    500 um length scale is present in this data regardless of a, g.
  * box half-width in x,y = 5a + g  (scales with a)

Writes island_cap_3d.png and island_cap_grid.csv
"""
import os, json, glob, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
EPS0 = 8.8541878128e-12
EPS_R = 11.45
eps_eff = EPS0 * (1 + EPS_R) / 2

rows = [json.load(open(f)) for f in sorted(glob.glob(f"{HERE}/coeff/results/*.json"))]
a = np.array([r["radius_um"] for r in rows])
g = np.array([r["gap_um"] for r in rows])
C = np.array([r["C_island_ground_fF"] for r in rows])

av = np.array(sorted(set(a))); gv = np.array(sorted(set(g)))
A, G = np.meshgrid(av, gv, indexing="ij")


def grid(vals):
    M = np.full_like(A, np.nan, float)
    for ai, gi, vi in zip(a, g, vals):
        M[np.where(av == ai)[0][0], np.where(gv == gi)[0][0]] = vi
    return M


Cg = grid(C)

fig = plt.figure(figsize=(20, 8.5))

specs = [
    (Cg,                 "C  (fF)",                     "viridis", (24, -128),
     "C(a, g)"),
    (Cg / A,             "C / a   (fF/µm)",             "plasma",  (24, -128),
     "C / a   (leading $C\\propto a$ divided out)"),
    (Cg / (EPS0 * A * 1e-6) / 1e15 * 1e15,  # = C[F]/(eps0 a[m]) ; keep dimensionless
     r"$C / (\varepsilon_0\, a)$", "cividis", (20, -110),
     r"$C/(\varepsilon_0 a)$   (dimensionless; isolated disc $\to 8\varepsilon_r{+}$...)"),
]
# fix panel-3 data properly
P3 = grid(C * 1e-15 / (EPS0 * a * 1e-6))

for i, (Z, zlab, cmap, view, title) in enumerate(specs, 1):
    ax = fig.add_subplot(2, 3, i, projection="3d")
    Zp = P3 if i == 3 else Z
    ax.plot_surface(A, G, Zp, cmap=cmap, alpha=.9, edgecolor="k", lw=.35, antialiased=True)
    zpts = (C * 1e-15 / (EPS0 * a * 1e-6)) if i == 3 else (C if i == 1 else C / a)
    ax.scatter(a, g, zpts, c="r", s=28, depthshade=False)
    ax.set_xlabel("radius $a$ (µm)", labelpad=8)
    ax.set_ylabel("gap $g$ (µm)", labelpad=8)
    ax.set_zlabel(zlab, labelpad=8)
    ax.set_title(title, fontsize=10)
    ax.view_init(elev=view[0], azim=view[1])

# bottom row: 2D slices for reading off the functional form
ax = fig.add_subplot(2, 3, 4)
for gi, c in zip(gv, plt.cm.viridis(np.linspace(0, .85, len(gv)))):
    m = g == gi
    o = np.argsort(a[m])
    s, b = np.polyfit(a[m], C[m], 1)
    xs = np.linspace(a[m].min(), a[m].max(), 20)
    ax.plot(xs, s * xs + b, "-", color=c, lw=1)
    ax.plot(a[m][o], C[m][o], "o", color=c,
            label=f"g={gi:.0f}: C={s:.3f}·a{b:+.0f}")
ax.set_xlabel("radius $a$ (µm)"); ax.set_ylabel("C (fF)")
ax.set_title("C vs a is LINEAR at fixed g  (≤1.4%)\nslope s(g), tiny intercept")
ax.grid(alpha=.3); ax.legend(fontsize=8)

ax = fig.add_subplot(2, 3, 5)
for ai, c in zip(av, plt.cm.plasma(np.linspace(0, .85, len(av)))):
    m = a == ai
    o = np.argsort(g[m])
    ax.plot(g[m][o], C[m][o], "o-", color=c, label=f"a={ai:.0f} µm")
ax.set_xlabel("gap $g$ (µm)"); ax.set_ylabel("C (fF)")
ax.set_xscale("log")
ax.set_title("C vs g   (per radius, log-x)"); ax.grid(alpha=.3, which="both"); ax.legend(fontsize=8)

ax = fig.add_subplot(2, 3, 6)
L = np.log(a / g)
Y = C * 1e-15 / (eps_eff * a * 1e-6)
for ai, c in zip(av, plt.cm.plasma(np.linspace(0, .85, len(av)))):
    m = a == ai
    o = np.argsort(L[m])
    ax.plot(L[m][o], Y[m][o], "o-", color=c, label=f"a={ai:.0f} µm")
ax.set_xlabel(r"$\ln(a/g)$"); ax.set_ylabel(r"$C/(\varepsilon_{\rm eff}\,a)$")
ax.set_title(r"dimensionless vs $\ln(a/g)$" "\n(curves offset by $a$ = broken $a/g$ scaling; concave)")
ax.grid(alpha=.3); ax.legend(fontsize=8)

fig.suptitle(f"Circular island → ground capacitance, 4×4 (a,g) grid   "
             f"($\\varepsilon_r$={EPS_R}, order-2, 500 µm grounded box top/bottom)", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(f"{HERE}/island_cap_3d.png", dpi=140)
print("wrote island_cap_3d.png")

with open(f"{HERE}/island_cap_grid.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["a_um", "g_um", "a_over_g", "C_fF", "C_per_a_fF_per_um",
                "C_over_eps0_a", "C_over_epseff_a", "ln_a_over_g"])
    for i in np.argsort(a / g):
        w.writerow([f"{a[i]:.0f}", f"{g[i]:.0f}", f"{a[i]/g[i]:.3f}", f"{C[i]:.3f}",
                    f"{C[i]/a[i]:.4f}", f"{C[i]*1e-15/(EPS0*a[i]*1e-6):.3f}",
                    f"{Y[i]:.3f}", f"{np.log(a[i]/g[i]):.3f}"])
print("wrote island_cap_grid.csv")
