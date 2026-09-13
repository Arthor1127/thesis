# `ansys_results/` — the eigenmode-Qc cross-check, handed back

Read **`SUMMARY.md`** first — it is the answer to `../HANDOFF.md`.

## What's here

```
SUMMARY.md                       the written answer  <- start here
qc_vs_R.png                      Qc vs R and ratio-to-analytic, all methods

bareres_R*_result.json           HFSS DRIVEN S11 fit (the primary deliverable):
                                   fr_ghz, Qc_oneport_fit, Qi_oneport_fit,
                                   Qc_group_delay, sweep window, convergence tail
bareres_R*.s1p                    raw driven S11 vs frequency (Touchstone,
                                   reference impedance = R, renormalize off)
bareres_R*_s11.csv               same data as CSV: freq_hz, S11_mag_db,
                                   S11_phase_deg, S11_re, S11_im
bareres_R*_convergence.csv        adaptive-pass history (pass, tets, deltaS)

bareres_eig_R*_result.json       HFSS EIGENMODE result per R (secondary):
                                   fr, complex eigenfrequency, Qc = Re(f)/(2 Im f)

run_bareres_s11.py               the driven-S11 driver
bareres_eigen.py                 the eigenmode driver
make_figure.py                   regenerates qc_vs_R.png from the *_result.json
_fit_offline.py                  re-runs the S11 fit from a saved .s1p (used for
                                   R=5e6, whose driver crashed post-export)
aedt_projects/                   the .aedt project files
```

`bareres_eigen.py` and `run_bareres_s11.py` import shared HFSS helpers from
`../../ansys_impedance_handoff/ansys_results/` (`run_ansys_tests.py`,
`gds_geometry.py`) — keep that sibling directory in place.

## Status per R — all three done

| R (Ω) | driven S11 `fr` / `Qc` | eigenmode `Qc` |
|---|---|---|
| 2e5 | 7.98931 GHz / 3790  | 2529  |
| 1e6 | 7.98931 GHz / 18904 | 11698 |
| 5e6 | 7.98931 GHz / 94457 | 57883 |

Driven S11 `Qc` is a flat 1.20 × analytic across all three R. It matches
PALACE at R ≤ 1e6; at R = 5e6 PALACE drops to ≈ analytic while HFSS driven
does not — see `SUMMARY.md`.

## Running

In the `qiskit_qc` conda env, from this directory:

```
python run_bareres_s11.py            # driven S11, all three R
python run_bareres_s11.py 1e6        # one R
python bareres_eigen.py              # eigenmode Qc, all three R
python make_figure.py                # rebuild qc_vs_R.png
```

Env knobs (driven): `HFSS_CORES` (default 6), `BR_MINPASS` / `BR_MAXPASS`
(forced adaptive passes — delta-S convergence is otherwise blind to a
weakly-loaded resonance), `BR_FS` (sweep centre frequency, GHz — must sit on
the resonance or the sweep never samples the pole), `BR_SPAN_LW` (half-window
in linewidths), `BR_SWEEP` (`Interpolating` or `Discrete`), `BR_NPTS`,
`BR_FADAPT` (off-resonance adapt frequency), `SI_TANDELTA` (default 1.2e-5).
Eigenmode: `BR_DF` (MaxDeltaFreq %), `SI_TANDELTA` forced to 0.

## Key setup choices

- The lumped element (driven port impedance, or eigenmode parallel-RLC R) sits
  on a sheet spanning **one** CPW gap, trace-to-ground. A sheet spanning both
  gaps puts two R's in parallel and halves the extracted Qc.
- Driven: adaptive mesh solved **off** resonance (uniform line → frequency-
  independent transverse mode), which dodges the near-singular FEM system at a
  high-Q resonance; a small `tan δ` gives `|S11|` a real, fittable dip; the
  complex one-port fit (Probst et al. 2015) then separates `Qc` from `Qi`.
- Eigenmode: lossless silicon + PEC metal, so `Q = Qc` exactly. Forced
  16–26 adaptive passes to `MaxDeltaFreq = 0.02 %`.
