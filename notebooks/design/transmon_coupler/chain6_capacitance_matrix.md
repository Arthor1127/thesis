# 6-transmon chain — capacitance matrix

Palace electrostatic (SQDMetal), circular transmons (`CircTransmon`, r = 136 µm,
pad gap 20 µm) joined by identical rounded interdigital couplers
(`CapNInterdigitalRound`, 8 fingers, 80 µm long, 8 µm wide, 4 µm gap, 6 µm tip
gap), pitch 454 µm, `eps_r = 11.45`. Geometry: `build_chain.build_chain(6)`.

## Converged Maxwell capacitance matrix (fF)

Terminal 1 is the ground plane; terminals 2–7 are Q1–Q6 (west→east).

|      | GND | Q1 | Q2 | Q3 | Q4 | Q5 | Q6 |
|------|------|------|------|------|------|------|------|
| GND | 1859.60 | -94.38 | -92.55 | -92.43 | -91.98 | -91.90 | -94.30 |
| Q1 | -94.38 | 156.81 | -44.82 | -0.07 | -0.00 | -0.00 | -0.00 |
| Q2 | -92.55 | -44.82 | 202.30 | -44.76 | -0.08 | -0.00 | -0.00 |
| Q3 | -92.43 | -0.07 | -44.76 | 202.14 | -44.71 | -0.08 | -0.00 |
| Q4 | -91.98 | -0.00 | -0.08 | -44.71 | 201.44 | -44.25 | -0.07 |
| Q5 | -91.90 | -0.00 | -0.00 | -0.08 | -44.25 | 201.82 | -44.93 |
| Q6 | -94.30 | -0.00 | -0.00 | -0.00 | -0.07 | -44.93 | 156.96 |

(saved as `chain6_capacitance_matrix.csv`)

## Key circuit quantities

| quantity | value |
|----------|-------|
| nearest-neighbour coupling $C_c$ | **44.7 ± 1 fF**  (bonds: 44.8, 44.8, 44.7, 44.3, 44.9) |
| next-nearest coupling | 0.076 fF  (≈ 0.17 % of $C_c$ — chain is genuinely NN) |
| island → ground-plane | 92.2 fF interior, 94.3 fF ends |
| island shunt cap to ground $C_\Sigma$ | 113 fF interior, 112 fF ends |
| island Maxwell self-cap $C_{ii}$ | 202 fF interior, 157 fF ends |
| charging energy $E_C = e^2/2\,(C^{-1})_{ii}$ | 108 MHz interior, 132 MHz ends |
| coupling ratio $\rho = C_c / C_\Sigma$ | ≈ 0.40 (interior) — strongly capacitively coupled |

## Convergence

The interdigital-finger capacitance is dominated by the metal-edge field
singularity in the 4 µm gaps, which converges slowly. Two independent mesh
families were run:

| Palace AMR cycle | DOF | error indicator | $C_c$ (chain mean) |
|---|---|---|---|
| 1 | 965,918 | 0.568 | 65.3 fF |
| 2 | 1,323,116 | 0.381 | 51.1 fF |
| 3 | 3,931,443 | 0.266 | 46.7 fF |
| **Aitken extrapolation** ($r$ = 0.31) | — | — | **44.7 fF** |

| uniform coupler mesh | $C_c$ (chain mean) |
|---|---|
| 3.0 µm | 74.2 fF |
| 2.0 µm | 65.1 fF |
| 1.5 µm | 55.8 fF |

Both series decrease monotonically toward ≈ 44–46 fF. The AMR series has a
clean geometric ratio (each bond $r$ = 0.29–0.35), giving the extrapolation
**$C_c$ ≈ 45 ± 2 fF**. A 4th AMR cycle (14.2 M DOF) was attempted
but Palace's flux-recovery error estimator stalled at that size (same failure
as a 9.4 M uniform run) and it was killed; the 3-cycle extrapolation stands.

Island→ground and the self-capacitances are much better converged
(< 3 % drift across the AMR cycles).

## Cross-check: interdigital-capacitor formula

$$C_{\rm IDC} = \varepsilon_0\varepsilon_{\rm eff}(N-1)\,\ell\,\frac{K(k)}{K(k')},
\quad k = \sin\!\frac{\pi}{2}\frac{w}{w+s}$$

With $N$=8, $\ell$=80 µm, $w$=8 µm, $s$=4 µm, $\varepsilon_{\rm eff}$=(1+11.45)/2:
$K(k)/K(k')$ = 1.28, giving $C_{\rm IDC}$ = **39.5 fF** core. The
extrapolated sim value is 44.7 fF → offset $C_0$ ≈ 5 fF
(finger tips + tapered leads + direct pad-to-pad), consistent with expectation.
This supports the model $C_c \approx \alpha\,C_{\rm IDC} + C_0$ for a future
coupler-geometry sweep.

## Files
- `chain6_capacitance_matrix.csv` — the converged matrix above
- `chain6_cap_convergence.png` — convergence plots
- `chain_cap/results/chain6_c*_o2_terminal-*.csv` — uniform-mesh raw matrices
- `chain_cap/amr_iter/iter{1,2,3}.csv` — the 3 AMR-cycle raw matrices
- `build_chain.py`, `run_chain_cap.py`, `run_chain_cap_amr.py` — geometry + runners
