# Transmon interdigital coupler — coupling-capacitance model

Two grounded `CircTransmon`s (r = 136 µm, pad gap 20 µm) joined by one
`CapNInterdigitalRound` coupler, each comb a galvanic extension of one island.
The pads are held **60 µm of bare lead apart** on each side so the measured
island–island mutual `|C_Maxwell[Q1][Q2]|` is the coupler capacitance `C_c`
(direct pad-to-pad is 0.3–0.8 fF, from no-coupler reference runs).
`eps_r = 11.45`. Geometry: `build_coupler_pair.build_coupler_pair()`.

Each point: Palace electrostatic + adaptive mesh refinement, `C_c` Aitken-
extrapolated from the per-cycle sequence (4–6 cycles; the ■ points reached only 3).

## Model

$$\boxed{\;C_c \;\approx\; \alpha\,C_{\rm IDC} \;+\; C_0\;}$$

$$C_{\rm IDC} = \varepsilon_0\,\varepsilon_{\rm eff}\,(N-1)\,\ell_{\rm eff}\;\frac{K(k)}{K(k')},
\qquad \ell_{\rm eff} = \ell - t,
\qquad k = \sin\!\Big(\frac{\pi}{2}\frac{w}{w+s}\Big),
\qquad \varepsilon_{\rm eff} = \frac{1+\varepsilon_r}{2}$$

N = finger count, ℓ = finger length, **ℓ_eff = ℓ − t = the finger overlap**
(t = finger_tip_gap), w = finger width, s = finger gap. K = complete elliptic
integral of the first kind.

**Fit (15 points, weighted by AMR-cycle count):**

| | α | C₀ (fF) | RMS resid | max resid |
|---|---|---------|-----------|-----------|
| **`α·C_IDC(ℓ−t) + C₀`** | **1.021** | **6.56** | **0.87 fF** | **1.39 fF** |
| + explicit tip term (β·C_tip) | 0.987 | 5.19 | 0.43 fF | 0.93 fF |

(the second row adds `C_tip = ε₀ε_eff·N·w·K(k_t)/K(k_t')`, k_t from t, β = 0.64;
it barely helps once ℓ_eff already carries the tip-gap dependence.)

**Reading the coefficients:**

- **α ≈ 1.02 — the conformal interdigital-capacitor formula is essentially
  exact for the fingers** (within ~2 %), *provided* the active length is the
  finger overlap ℓ − t, not the full finger. Coupling comes from where the two
  combs' fingers run alongside each other; the tip-to-busbar stub barely
  contributes. Using the full ℓ instead gives α ≈ 0.98, C₀ ≈ 5 fF but RMS
  1.2 fF and a structured tip-gap residual up to 3.5 fF.
- **C₀ ≈ 7 fF — a geometry-fixed offset**: the finger tips facing the
  opposing bus bar, the two tapered necks running antiparallel from the combs
  to the pads, and ~0.5 fF direct pad-to-pad. Roughly constant across the
  sweep because the neck length is held fixed (60 µm).
- Validity: fitted over N = 4–16, ℓ = 40–160 µm, s = 3–10 µm, w = 6–12 µm,
  t = 3–12 µm → C_c = 21–88 fF. Cross-check: the same coupler (N8/L80/g4/w8/t6)
  in the 6-transmon chain gave C_c = 44.7 fF; here 43.8 fF; model 43.8 fF.

## Sweep data (Aitken-extrapolated C_c)

| N | ℓ | s | w | t | C_IDC (fF) | C_c (fF) | AMR cyc |
|---|---|---|---|---|-----------|----------|---------|
| 4 | 80 | 4 | 8 | 6 | 15.7 | 21.2 | 4 |
| 8 | 40 | 4 | 8 | 6 | 16.8 | 25.1 | 4 |
| 6 | 80 | 4 | 8 | 6 | 26.1 | 32.4 | 4 |
| 8 | 80 | 10 | 8 | 6 | 26.4 | 34.1 | 4 |
| 8 | 80 | 6 | 8 | 6 | 31.7 | 39.4 | 4 |
| 8 | 80 | 4 | 6 | 6 | 33.0 | 39.6 | 4 |
| 8 | 80 | 4 | 8 | 12 | 33.6 | 40.3 | 4 |
| 8 | 80 | 4 | 8 | 6 | 36.5 | 43.8 | 4 |
| 8 | 80 | 4 | 8 | 3 | 38.0 | 46.1 | 5 |
| 8 | 80 | 3 | 8 | 6 | 40.3 | 47.5 | 4 |
| 8 | 80 | 4 | 12 | 6 | 41.9 | 50.0 | 5 |
| 8 | 120 | 4 | 8 | 6 | 56.3 | 63.0 | 5 |
| 12 | 80 | 4 | 8 | 6 | 57.4 | 65.9 | 5 |
| 8 | 160 | 4 | 8 | 6 | 76.0 | 82.9 | 6 |
| 16 | 80 | 4 | 8 | 6 | 78.3 | 87.6 | 5 |

## Files
- `fit_coupler_model.py` — this fit
- `coupler_model_fit.png` — C_c vs each swept parameter, with the model
- `coupler_model_fit.json` — coefficients + per-point data
- `build_coupler_pair.py`, `run_coupler_pair_cap.py` — geometry + AMR runner
- `coupler_sweep/results/*_amr.json` — per-point raw matrices + AMR cycles
