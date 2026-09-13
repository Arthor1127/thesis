# ANSYS HFSS driven-S11 cross-check of PALACE's eigenmode-Qc method

Request: `../HANDOFF.md`. Geometry: a straight shorted quarter-wave CPW
resonator (`../geometry/bare_resonator.gds`, w = 20 µm, gap = 11.5 µm, length
3.773 mm, on 500 µm silicon εr = 11.45), PEC short at +x (STG), a **lumped port
of resistance R across one open-end CPW gap** at −x (OTG1). Analytic reference
`Qc = (π/4)(R/Z0)`, Z0 = 49.961 Ω.

## Bottom line

The HFSS **driven S11** fit gives, for all three R values, a coupling quality
factor that is a **flat 1.20 × the analytic `(π/4)(R/Z0)`**:

| R (Ω) | analytic `Qc` | PALACE eigenmode | **HFSS driven S11** | HFSS ÷ analytic | HFSS ÷ PALACE |
|---|---|---|---|---|---|
| 2e5 | 3144  | 3841  (×1.22) | **3790**  | 1.21 | 0.99 |
| 1e6 | 15720 | 18495 (×1.18) | **18904** | 1.20 | 1.02 |
| 5e6 | 78601 | 78182 (×0.99) | **94457** | 1.20 | 1.21 |

`fr` = **7.98931 GHz** for all three R (a real parallel R adds no reactance),
reproduced to 6 significant figures across the three independent driven runs.

**Against the handoff's two scenarios:**

- **At R = 2e5 and R = 1e6** the HFSS driven-S11 fit lands on **PALACE's
  numbers, not the analytic column** (0.99× and 1.02× PALACE; 1.21× and 1.20×
  analytic). This is the "reassuring cross-validation" outcome — two fully
  independent tools *and* techniques (HFSS driven S-parameter fit vs PALACE
  complex-eigenfrequency) agree to ≈ 2 %, and both say the analytic formula is
  ~18 % low here.

- **At R = 5e6 they diverge.** HFSS driven S11 stays on its flat 1.20× trend
  (Qc ≈ 94 500); PALACE eigenmode drops to 0.99× analytic (Qc ≈ 78 200). By
  the handoff's own logic — PALACE matching the analytic column at R = 5e6 but
  not at lower R — **PALACE's eigenmode method appears to carry an R-dependent
  bias that shows up at high Qc**. R = 5e6 is exactly the hardest case for an
  eigenmode solver: the mode's fractional damping is `Im(f)/Re(f) ≈ 8.6e-6`, so
  Qc comes from a part-per-100 000 imaginary component of the complex
  eigenfrequency. (HFSS's *own* eigenmode solve shows a different but also
  R-dependent offset — see the secondary table below — which is consistent with
  eigenmode Q-extraction being the fragile step here, independent of tool.)

The self-consistency check on the driven fit: it returns `Qi` = **92 200 –
92 500 in all three fits** — the `1/(q·tan δ)` dielectric limit for the
`tan δ = 1.2e-5` put on the silicon — even though `Qi` ranges from 25× to 1×
of `Qc` across the set. The complex one-port fit is cleanly separating the
port (coupling) loss from the dielectric loss. At R = 5e6 `Qi ≈ Qc`
(near-critical), where that separation is least well-conditioned; but pinning
`Qi = 92 400` from the two well-separated fits and back-solving `Qc` from the
measured `Ql = 46 659` gives `Qc = 94 260` — the same answer, so the ≈ 94 500
figure does not depend on the fit resolving the near-degeneracy.

## The S11 data

| R (Ω) | notch depth | `Ql` (fit) | `Qc` (fit) | `Qi` (fit) | sweep type |
|---|---|---|---|---|---|
| 2e5 | −0.7 dB  | 3641  | 3790  | 92 518 | interpolating, ±0.36 GHz |
| 1e6 | −3.6 dB  | 15693 | 18904 | 92 376 | interpolating, ±71 MHz |
| 5e6 | −38.4 dB | 46659 | 94457 | 92 206 | interpolating, ±2.5 MHz, centred on `fr` |

`Qc` / `Qi` / `Ql` from a complex one-port (reflection) resonance fit — Probst
et al. 2015: `S11 = env·[1 − (2Ql/|Qc|)e^{−iφ}/(1 + 2iQl(f/fr−1))]`. The
group-delay estimate `Qc = π fr τg/2` agrees for R = 2e5 / 1e6 (3781, 19694)
but is unusable at R = 5e6 where the −38 dB (near-critical) notch makes the
phase non-monotonic — the one-port fit is the reliable estimator there.

R = 5e6 needed care: a fine adaptive mesh (16–20 passes, MaxΔS = 1e-4) *and*
an interpolating sweep window only ±25 linewidths wide, centred exactly on the
already-known `fr`, so the rational-model sweep's first samples straddle the
0.1 MHz-wide pole. Wider windows or a coarser mesh return a flat `S11`; an
every-point discrete sweep dense enough to resolve the pole exhausts memory on
this machine (~5 GB free). The `.s1p` was solved cleanly; only the driver's
post-processing crashed, so its fit was re-run offline (`_fit_offline.py`).

## HFSS eigenmode (secondary — cruder port model)

An HFSS eigenmode solve, lossless dielectric, lumped **parallel-RLC (R only)**
on the same one-gap sheet:

| R (Ω) | HFSS eigenmode `Qc` | ×analytic | `fr` |
|---|---|---|---|
| 2e5 | 2529  | 0.80 | 7.976 GHz |
| 1e6 | 11698 | 0.74 | 7.978 GHz |
| 5e6 | 57883 | 0.74 | 7.978 GHz |

Consistently **~0.75× analytic**, and `fr` ~13 MHz below the driven value —
i.e. the lumped-RLC element on a 10 µm partial-gap sheet is not the same
boundary as a real port across the full CPW gap. Treated as a rough
consistency check on the power law (`Qc ∝ R`, exponent 0.97), not as an
absolute number; the driven S11 result is the one to trust.

## Files

- `bareres_R{2e05,1e06,5e06}_result.json` — driven fit (the deliverable):
  `fr_ghz`, `Qc_oneport_fit`, `Qi_oneport_fit`, `Ql_oneport_fit`,
  `Qc_group_delay`, sweep window (+ adaptive-convergence tail for R = 2e5/1e6).
- `bareres_R*.s1p`, `bareres_R*_s11.csv` — raw driven S11 vs frequency
  (Touchstone reference impedance = R, `renormalize = False`).
- `bareres_R{2e05,1e06}_convergence.csv` — adaptive-pass history (the R = 5e6
  driver crashed before exporting its own; same setup as R = 1e6, which
  converged to ΔS = 6e-5 by pass 16).
- `bareres_eig_R*_result.json` — HFSS eigenmode (secondary).
- `qc_vs_R.png` — Qc vs R and ratio-to-analytic, all four estimates.
- `run_bareres_s11.py`, `bareres_eigen.py`, `make_figure.py`, `_fit_offline.py`
  — drivers.
- `aedt_projects/` — the `.aedt` project files.

## Method detail (HFSS)

ANSYS Electronics Desktop 2024 R2 via `ansys-aedt-core` 1.5 (COM transport).
Geometry: clean rectangles from the GDS trace bbox; PEC ground with the CPW
slot, grounded PEC backplane 500 µm below the trace, compact vacuum box,
zero-thickness Perfect-E metal. Driven: lumped port `impedance = R` on a
one-gap sheet at the open end (a sheet across *both* gaps loads the mode with
R/2 and halves `Qc`); `renormalize = False` so the S11 reference is R itself;
a small `tan δ = 1.2e-5` on the silicon so `|S11|` has a real, fittable dip
rather than a pure phase all-pass; adaptive mesh solved **off** resonance at
8.30 GHz (for a uniform CPW line the transverse mode is frequency-independent),
which dodges the near-singular FEM system at a high-Q resonance; interpolating
(rational-model) sweep across the resonance.
