# Circular transmon island → ground: empirical capacitance model

## Geometry

A single circular metal disc of radius `a` centred in an infinite ground
plane, separated from it by an annular vacuum gap of width `g` (ground-plane
hole radius `b = a + g`). No junction, no coupler — just the island-to-ground
self-capacitance `C`.

```
        ┌───────────────── ground plane ─────────────────┐
        │            ╭───────────────╮                    │
        │          ╱   ╭─────────╮     ╲   ← gap g         │
        │         │   │  island   │     │                  │
        │         │   │  radius a │     │                  │
        │          ╲   ╰─────────╯     ╱                   │
        │            ╰───────────────╯                     │
        └────────────────────────────────────────────────┘
```

## Simulation setup (SQDMetal + Palace, electrostatic)

| item | value |
|------|-------|
| solver | Palace electrostatic, **order 2**, BoomerAMG/CG, tol 1e-8 |
| substrate permittivity | **ε_r = 11.45** (SQDMetal's "silicon"; *not* 11.9) |
| metal | zero-thickness sheet (quasi-3D); C is thickness-independent to <0.1 % over t = 30–800 nm |
| mesh | thin fine ring at the rim, `res_min = max(g/6, 1.2) µm` (~6 elements across the gap); mesh-converged (cross-check stable to 0.3 % between res 2 µm and 3 µm) |
| domain box | half-width `5a + g` in x,y; **±500 µm in z**, and **all six box faces are grounded PEC** (see caveats) |
| data | job 1938614, 16 points: `a ∈ {60, 100, 150, 220} µm` × `g ∈ {6, 12, 24, 45} µm`, ParaView output disabled |

## The model

The simulated `C(a, g)` is, to high accuracy, **linear in the radius at
fixed gap**:

$$\boxed{\;C(a,g) \;=\; s(g)\,a \;+\; b(g)\;}$$

with

$$s(g) = s_0 + s_1\ln g \;=\; 1.472 - 0.2235\,\ln g
\qquad [\text{fF}/\mu\text{m}]$$

$$b(g) = b_0 + b_1\sqrt{g} \;=\; -19.21 + 1.882\,\sqrt{g}
\qquad [\text{fF}]$$

**`a`, `g` in micrometres, `C` in femtofarads, `ln` natural, for ε_r = 11.45.**

Equivalent single line:

```
C_fF = (1.472 - 0.2235*log(g_um)) * a_um  +  (-19.21 + 1.882*sqrt(g_um))
```

Fitted jointly to all 16 points by least squares:

| coeff | value | 1σ |
|-------|-------|-----|
| s₀ | +1.4717 fF/µm | 0.012 |
| s₁ | −0.2235 fF/µm | 0.004 |
| b₀ | −19.21 fF | 1.32 |
| b₁ | +1.882 fF/µm^(1/2) | 0.28 |

**Accuracy over the fitted range: 2.4 % max, 0.9 % rms.**

### Interpretation of the terms

* **`s(g)·a`** — the dominant term. `dC/da` is set by the field near the
  island rim. It grows as the gap narrows because a tighter coplanar ground
  captures more field lines per unit rim length. Writing
  `s(g) = 0.2235·ln(725 µm / g)`, an apparent cutoff length ~700 µm appears —
  of order the box height, i.e. the coplanar logarithm is cut off by the
  finite grounded box rather than running forever (see caveats).
  For reference the *isolated* disc `C = 8ε_eff·a` has slope `8ε_eff ≈ 0.44
  fF/µm`; every `s(g)` here is 1.4–2.4× larger — the coplanar ground roughly
  doubles-to-triples `dC/da`.
* **`b(g)`** — a small negative intercept, −14 fF (g=6) rising to −7 fF
  (g=45). `C(a)` at fixed `g` has slight positive curvature, so its large-`a`
  tangent has a negative intercept; a wider gap weakens the edge effect and
  pushes the behaviour back toward the pure `C ∝ a` disc law, shrinking
  `|b|`. This is the least-certain part of the model: `b` is an
  extrapolation of a straight line fitted over `a = 60–220 µm` down to
  `a = 0`, and its 1σ bar is ±0.4–1 fF. `√g`, `ln g` and linear-in-`g` fits
  are statistically indistinguishable with 4 gaps; `√g` tracks the points
  best (chosen here).

## Validity

* Fitted range: `a = 60–220 µm`, `g = 6–45 µm`, i.e. `a/g = 1.3–37`.
  Interpolation inside this box is good to ~2 %.
* Cross-check: at `a = 135.94 µm, g = 20 µm` (the radius the earlier
  `8ε_eff·a·[1+(1/π)ln(a/g)]` estimate predicted for 100 fF) this model
  returns `C_model = 98.2 fF`; independent Palace runs gave 96.9–97.2 fF
  (~1 % — mesh-converged). Consistent.
* Do **not** extrapolate `b(g)` past the range; `b → 0` around `g ≈ 100 µm`
  in this parametrisation, which is not physical, just the end of the fit.

## Why not the textbook forms

| model | form | max error, all 16 pts |
|-------|------|----------------------|
| coaxial edge | `C = 2πε_eff·t / ln(b/a)` | 100–800× too small (it is only the metal edge-to-edge term; a <1 % correction) |
| isolated disc | `C = 8ε_eff·a` | −45 % (ignores the coplanar ground entirely) |
| disc + log gap | `C = 8ε_eff·a·[1 + (1/π)ln(a/g)]` | 15 % |
| best `A + B ln(a/g)` fit | `C/(ε_eff a) = 7.4 + 2.83 ln(a/g)` | 14 % |

**`C/(ε_eff·a)` is not a function of `a/g` alone.** Regressing
`C/a = α + β ln(a/g) + γ ln a` leaves `γ = −0.09`, half the size of `β` and
clearly non-zero: there is a genuine absolute-size dependence that the
`ln(a/g)` scaling cannot absorb, and the `C = s(g)a + b(g)` decomposition
captures it directly and far more accurately.

## Caveats

1. **ε_r = 11.45** is hard-wired into the coefficients. For a different
   dielectric, the leading behaviour scales roughly as `ε_eff = ε₀(1+ε_r)/2`,
   but this was not verified with a sweep — re-fit if ε_r changes.
2. **The Palace domain box has all six faces grounded**, with the top and
   bottom at ±500 µm from the island plane. This introduces a fixed ~500 µm
   length scale that is *not* part of an infinite-substrate coplanar problem
   — it is very likely why `C/(ε_eff a)` fails to collapse onto `a/g` (the
   `γ ln a` residual above) and why the `s(g)` logarithm appears to saturate
   at ~700 µm. For a real chip with a much thicker or open backside the
   size dependence would differ. A dedicated box-size / substrate-thickness
   sweep is needed to separate physical from boundary effects. Until then,
   treat this model as *specific to this simulation domain*, accurate for
   islands well inside the fitted `(a, g)` box.
3. `b(g)` is an extrapolated intercept (see above) — the single largest
   source of model uncertainty at small `a`.

## In code

Embedded in `parametric_coupling.DesignHelper`:

| method | does |
|--------|------|
| `circular_island_capacitance(radius, gap, eps_r=11.45)` | forward model → C |
| `circular_island_radius_for_capacitance(C, gap, eps_r=11.45)` | inverse for `a` (closed form) |
| `circular_island_gap_for_capacitance(C, radius, eps_r=11.45)` | inverse for `g` (bisection) |
| `circular_pad_capacitance(...)` | the old coaxial-edge formula — kept for reference only, do not use for a thin-film pad |

## Files

* `run_cap_island.py` — the per-point Palace runner
* `coeff/results/*.json` — the 16 raw results (job 1938614)
* `island_cap_grid.csv` — tabulated `C`, `C/a`, `C/(ε₀a)`, `ln(a/g)`
* `plot3d_coeff.py` → `island_cap_3d.png` — 3D surfaces + 2D slices
* `plot_bg.py` → `island_cap_bg.png` — `s(g)`, `b(g)` with error bars
* `fit_coeff.py` → `island_cap_coeff.{csv,png,json}` — the `ln(a/g)` fits (for the record)
