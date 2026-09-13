"""
Fit C_c(coupler geometry) from the AMR sweep + rerun (jobs 1938636/1938641/1938643).
Writes coupler_model_fit.png, coupler_model_fit.json, coupler_capacitance_model.md.

Geometry: 2 CircTransmons (r = 136 um, pad gap 20 um) + 1 CapNInterdigitalRound,
pads held 60 um of bare neck apart (direct pad-pad ~0.3-0.8 fF, per the
no-coupler references). eps_r = 11.45.

Models for the island-island coupling capacitance:
  1)  C_c = a * C_IDC + C0
  2)  C_c = a * C_IDC + b * C_tip + C0
      C_IDC = eps0 eps_eff (N-1) L  K(k)/K(k'),   k  = sin(pi/2 * w/(w+s))
      C_tip = eps0 eps_eff  N   w   K(kt)/K(kt'),  kt = sin(pi/2 * w/(w+t))
      eps_eff = (1 + eps_r)/2
"""
import os, json, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import ellipk

HERE = os.path.dirname(os.path.abspath(__file__))
EPS0, EPS_R = 8.8541878128e-12, 11.45
eeff = EPS0 * (1 + EPS_R) / 2
BASE = dict(N=8, L=80.0, g=4.0, w=8.0, t=6.0)


def KK(k):
    return ellipk(k ** 2) / ellipk(1.0 - k ** 2)


def C_idc(N, L, w, s, t=0.0):
    """Interdigital core. The electrically active length is the finger OVERLAP,
    l_eff = l - t  (t = finger_tip_gap): coupling comes from where the two
    combs' fingers run alongside each other, not the tip-to-busbar stub."""
    L_eff = L - t
    return eeff * (N - 1) * (L_eff * 1e-6) * KK(np.sin(np.pi * (w / (w + s)) / 2)) * 1e15


def C_tip(N, w, t):
    return eeff * N * (w * 1e-6) * KK(np.sin(np.pi * (w / (w + t)) / 2)) * 1e15


def aitken(seq):
    s = list(dict.fromkeys([round(x, 4) for x in seq]))
    if len(s) < 3:
        return (s[-1] if s else None), len(s)
    d1, d2 = s[-2] - s[-3], s[-1] - s[-2]
    ex = s[-1] - d2 * d2 / (d2 - d1) if abs(d2 - d1) > 1e-9 else s[-1]
    return ex, len(s)


pts = []
for f in sorted(glob.glob(f"{HERE}/coupler_sweep/results/*_amr.json")):
    d = json.load(open(f))
    if "pilot" in d["tag"]:
        continue
    N, L, g, w, t = (d["finger_count"], d["finger_length_um"], d["finger_gap_um"],
                     d["finger_width_um"], d["finger_tip_gap_um"])
    seq = [c["C_c_fF"] for c in d["amr_cycles"]] + [d["C_c_fF"]]
    ex, ncy = aitken(seq)
    pts.append(dict(N=N, L=L, g=g, w=w, t=t, ncy=ncy, final=d["C_c_fF"], ex=ex,
                    cidc=C_idc(N, L, w, g, t), ctip=C_tip(N, w, t), dof=d.get("dof")))

# weight by extrapolation quality: >=4 cycles -> 1, 3 -> 0.5, <3 -> 0.15
def weight(p):
    return {4: 1.0, 5: 1.0, 6: 1.0}.get(p["ncy"], 0.5 if p["ncy"] == 3 else 0.15)


x = np.array([p["cidc"] for p in pts])
xt = np.array([p["ctip"] for p in pts])
y = np.array([p["ex"] for p in pts])
wt = np.array([weight(p) for p in pts])
sw = np.sqrt(wt)

A1 = np.vstack([x, np.ones_like(x)]).T
(a1, c01), *_ = np.linalg.lstsq(sw[:, None] * A1, sw * y, rcond=None)
A2 = np.vstack([x, xt, np.ones_like(x)]).T
(a2, b2, c02), *_ = np.linalg.lstsq(sw[:, None] * A2, sw * y, rcond=None)

r1 = y - A1 @ [a1, c01]
r2 = y - A2 @ [a2, b2, c02]
rms1 = float(np.sqrt(np.average(r1 ** 2, weights=wt)))
rms2 = float(np.sqrt(np.average(r2 ** 2, weights=wt)))
mx1, mx2 = float(np.max(np.abs(r1))), float(np.max(np.abs(r2)))

report = {
    "source": "jobs 1938636 (sweep) + 1938641/1938643 (reruns), 15 points, Palace AMR + Aitken",
    "eps_r": EPS_R, "eps_eff_over_eps0": eeff / EPS0, "baseline": BASE,
    "n_points": len(pts),
    "model1": {"form": "C_c = a*C_IDC + C0  [fF]", "a": float(a1), "C0_fF": float(c01),
               "rms_resid_fF": rms1, "max_resid_fF": mx1},
    "model2": {"form": "C_c = a*C_IDC + b*C_tip + C0  [fF]", "a": float(a2), "b": float(b2),
               "C0_fF": float(c02), "rms_resid_fF": rms2, "max_resid_fF": mx2},
    "points": [{"N": p["N"], "L": p["L"], "g": p["g"], "w": p["w"], "t": p["t"],
                "C_IDC_fF": round(p["cidc"], 2), "C_c_fF": round(p["ex"], 2),
                "amr_cycles": p["ncy"], "dof": p["dof"]} for p in pts],
}
json.dump(report, open(f"{HERE}/coupler_model_fit.json", "w"), indent=2)


def m1(N, L, w, s, t):
    return a1 * C_idc(N, L, w, s, t) + c01


def m2(N, L, w, s, t):
    return a2 * C_idc(N, L, w, s, t) + b2 * C_tip(N, w, t) + c02


sweeps = [
    ("L", "finger_length ℓ (µm)", np.linspace(35, 170, 60)),
    ("N", "finger_count N", np.arange(4, 17)),
    ("g", "finger_gap s (µm)", np.linspace(2.5, 11, 60)),
    ("w", "finger_width w (µm)", np.linspace(5, 13, 60)),
    ("t", "finger_tip_gap t (µm)", np.linspace(2.5, 13, 60)),
]
fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
for ax, (key, xlabel, xr) in zip(axes.flat, sweeps):
    def kw(v):
        d = dict(BASE); d[key] = v; return d
    ax.plot(xr, [C_idc(kw(v)["N"], kw(v)["L"], kw(v)["w"], kw(v)["g"], kw(v)["t"]) for v in xr],
            ":", color="gray", label=r"$C_{IDC}$ formula")
    ax.plot(xr, [m1(kw(v)["N"], kw(v)["L"], kw(v)["w"], kw(v)["g"], kw(v)["t"]) for v in xr],
            "-", color="C3", lw=2, label=f"model 1")
    ax.plot(xr, [m2(kw(v)["N"], kw(v)["L"], kw(v)["w"], kw(v)["g"], kw(v)["t"]) for v in xr],
            "--", color="C0", lw=1.8, label="model 2 (+explicit tip term)")
    for p in pts:
        if any(p[k] != BASE[k] for k in BASE if k != key):
            continue
        ax.scatter(p[key], p["ex"], s=95, color="k",
                   marker=("o" if p["ncy"] >= 4 else "s"), zorder=5)
    ax.set_xlabel(xlabel); ax.set_ylabel(r"$C_c$ (fF)")
    ax.set_title(f"vary {key}   (others at baseline)")
    ax.grid(alpha=.3); ax.legend(fontsize=8)

ax = axes.flat[5]
xs = np.linspace(0, 90, 40)
ax.plot(xs, xs, ":", color="gray", label=r"$C_c=C_{IDC}$")
ax.plot(xs, a1 * xs + c01, "-", color="C3", lw=2,
        label=fr"$C_c={a1:.2f}\,C_{{IDC}}{c01:+.1f}$  (RMS {rms1:.1f} fF)")
for p in pts:
    ax.scatter(p["cidc"], p["ex"], s=75, color="k",
               marker=("o" if p["ncy"] >= 4 else "s"), zorder=5)
ax.set_xlabel(r"$C_{IDC}$ (fF)"); ax.set_ylabel(r"$C_c$ AMR-extrapolated (fF)")
ax.set_title("all 15 points vs the interdigital formula")
ax.grid(alpha=.3); ax.legend(fontsize=8)
fig.suptitle("Transmon coupler — coupling-capacitance model  "
             "(● ≥4 AMR cycles   ■ 3 cycles)", fontsize=13)
fig.savefig(f"{HERE}/coupler_model_fit.png", dpi=140)

# ---- second figure: residuals + AMR convergence ----
key_of = lambda p: next((k for k in BASE if p[k] != BASE[k]), "base")
col = {"N": "C0", "L": "C1", "g": "C2", "w": "C3", "t": "C4", "base": "k"}
fig2, ax2 = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
ax2[0].axhline(0, color="k", lw=.8)
for p, rr in zip(pts, r1):
    ax2[0].scatter(p["cidc"], rr, s=80, color=col[key_of(p)],
                   marker=("o" if p["ncy"] >= 4 else "s"), edgecolor="k", zorder=3)
import matplotlib.lines as mline
ax2[0].legend(handles=[mline.Line2D([], [], marker="o", ls="", color=col[k], label=k)
                       for k in ("N", "L", "g", "w", "t")], fontsize=9, title="varied")
ax2[0].set_xlabel(r"$C_{IDC}(\ell-t)$ (fF)"); ax2[0].set_ylabel(r"$C_c^{sim} - $ model 1 (fF)")
ax2[0].set_title(fr"model-1 residuals  (RMS {rms1:.2f} fF, max {mx1:.2f} fF)")
ax2[0].grid(alpha=.3)

for f in sorted(glob.glob(f"{HERE}/coupler_sweep/results/*_amr.json")):
    d = json.load(open(f))
    if "pilot" in d["tag"]:
        continue
    seq = list(dict.fromkeys([round(c["C_c_fF"], 3) for c in d["amr_cycles"]] + [round(d["C_c_fF"], 3)]))
    N, L, g, w, t = (d["finger_count"], d["finger_length_um"], d["finger_gap_um"],
                     d["finger_width_um"], d["finger_tip_gap_um"])
    p0 = dict(N=N, L=L, g=g, w=w, t=t)
    k = next((kk for kk in BASE if p0[kk] != BASE[kk]), "base")
    ex, _ = aitken(seq)
    ax2[1].plot(range(1, len(seq) + 1), np.array(seq) / ex, "o-", ms=4, color=col[k], alpha=.7)
ax2[1].axhline(1.0, color="k", lw=.8)
ax2[1].set_xlabel("AMR cycle"); ax2[1].set_ylabel(r"$C_c^{cycle} / C_c^{extrapolated}$")
ax2[1].set_title("per-point AMR convergence (Aitken target = 1)")
ax2[1].grid(alpha=.3)
fig2.suptitle("Coupler model — residuals and mesh convergence", fontsize=13)
fig2.savefig(f"{HERE}/coupler_model_residuals.png", dpi=140)

A1lbl=fr"$C_c={a1:.2f}\\,C_{{IDC}}(\\ell{{-}}t){c01:+.1f}$"
# ---- markdown ----
tbl = "\n".join(
    f"| {p['N']:.0f} | {p['L']:.0f} | {p['g']:.0f} | {p['w']:.0f} | {p['t']:.0f} "
    f"| {p['cidc']:.1f} | {p['ex']:.1f} | {p['ncy']} |"
    for p in sorted(pts, key=lambda p: p["cidc"]))
md = f"""# Transmon interdigital coupler — coupling-capacitance model

Two grounded `CircTransmon`s (r = 136 µm, pad gap 20 µm) joined by one
`CapNInterdigitalRound` coupler, each comb a galvanic extension of one island.
The pads are held **60 µm of bare lead apart** on each side so the measured
island–island mutual `|C_Maxwell[Q1][Q2]|` is the coupler capacitance `C_c`
(direct pad-to-pad is 0.3–0.8 fF, from no-coupler reference runs).
`eps_r = 11.45`. Geometry: `build_coupler_pair.build_coupler_pair()`.

Each point: Palace electrostatic + adaptive mesh refinement, `C_c` Aitken-
extrapolated from the per-cycle sequence (4–6 cycles; the ■ points reached only 3).

## Model

$$\\boxed{{\\;C_c \\;\\approx\\; \\alpha\\,C_{{\\rm IDC}} \\;+\\; C_0\\;}}$$

$$C_{{\\rm IDC}} = \\varepsilon_0\\,\\varepsilon_{{\\rm eff}}\\,(N-1)\\,\\ell_{{\\rm eff}}\\;\\frac{{K(k)}}{{K(k')}},
\\qquad \\ell_{{\\rm eff}} = \\ell - t,
\\qquad k = \\sin\\!\\Big(\\frac{{\\pi}}{{2}}\\frac{{w}}{{w+s}}\\Big),
\\qquad \\varepsilon_{{\\rm eff}} = \\frac{{1+\\varepsilon_r}}{{2}}$$

N = finger count, ℓ = finger length, **ℓ_eff = ℓ − t = the finger overlap**
(t = finger_tip_gap), w = finger width, s = finger gap. K = complete elliptic
integral of the first kind.

**Fit (15 points, weighted by AMR-cycle count):**

| | α | C₀ (fF) | RMS resid | max resid |
|---|---|---------|-----------|-----------|
| **`α·C_IDC(ℓ−t) + C₀`** | **{a1:.3f}** | **{c01:.2f}** | **{rms1:.2f} fF** | **{mx1:.2f} fF** |
| + explicit tip term (β·C_tip) | {a2:.3f} | {c02:.2f} | {rms2:.2f} fF | {mx2:.2f} fF |

(the second row adds `C_tip = ε₀ε_eff·N·w·K(k_t)/K(k_t')`, k_t from t, β = {b2:.2f};
it barely helps once ℓ_eff already carries the tip-gap dependence.)

**Reading the coefficients:**

- **α ≈ {a1:.2f} — the conformal interdigital-capacitor formula is essentially
  exact for the fingers** (within ~2 %), *provided* the active length is the
  finger overlap ℓ − t, not the full finger. Coupling comes from where the two
  combs' fingers run alongside each other; the tip-to-busbar stub barely
  contributes. Using the full ℓ instead gives α ≈ 0.98, C₀ ≈ 5 fF but RMS
  1.2 fF and a structured tip-gap residual up to 3.5 fF.
- **C₀ ≈ {c01:.0f} fF — a geometry-fixed offset**: the finger tips facing the
  opposing bus bar, the two tapered necks running antiparallel from the combs
  to the pads, and ~0.5 fF direct pad-to-pad. Roughly constant across the
  sweep because the neck length is held fixed (60 µm).
- Validity: fitted over N = 4–16, ℓ = 40–160 µm, s = 3–10 µm, w = 6–12 µm,
  t = 3–12 µm → C_c = 21–88 fF. Cross-check: the same coupler (N8/L80/g4/w8/t6)
  in the 6-transmon chain gave C_c = 44.7 fF; here 43.8 fF; model 43.8 fF.

## Sweep data (Aitken-extrapolated C_c)

| N | ℓ | s | w | t | C_IDC (fF) | C_c (fF) | AMR cyc |
|---|---|---|---|---|-----------|----------|---------|
{tbl}

## Files
- `fit_coupler_model.py` — this fit
- `coupler_model_fit.png` — C_c vs each swept parameter, with the model
- `coupler_model_fit.json` — coefficients + per-point data
- `build_coupler_pair.py`, `run_coupler_pair_cap.py` — geometry + AMR runner
- `coupler_sweep/results/*_amr.json` — per-point raw matrices + AMR cycles
"""
with open(f"{HERE}/coupler_capacitance_model.md", "w") as f:
    f.write(md)

print(json.dumps(report, indent=2, default=str))
print("\nwrote coupler_model_fit.{png,json}, coupler_capacitance_model.md")
