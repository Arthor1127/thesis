# `ansys_results/` — the ANSYS cross-check, handed back

Start with **`SUMMARY.md`** — the written answer to `../HANDOFF.md`.

Everything is in this one directory:

```
SUMMARY.md              the written answer  <- read this first
README.md               this file
comparison.png          all methods on one bar chart
testA_line_only.png     Test A: Zpi / Zin / S-params vs frequency
testC_full_design.png   Test C: driven S21 (flat) + eigenmode f_r
testC_coupler_geometry.png   why the hanger couples weakly
test{A,B,C}_*            raw data, one set per test (see table below)
testC_eigenmodes.json   the eigenmode result (resonator f_r)
*.py                    the scripts that produced all of it
aedt_projects/          the four .aedt project files (geometry + setup;
                        the large solved-field folders are left in
                        ~/Documents/Ansoft/ansys_test*.aedtresults/)
```

## Data, per test

`<tag>` is `testA_line_only`, `testA_yh150` / `_yh400` / `_yh1000` (ground-width
sweep), `testB_launchpads`, `testC_full_design`.

| file | contents |
|---|---|
| `<tag>.s2p` | raw HFSS Touchstone — **not renormalised**, reference = each port's field `Zpi(f)` |
| `<tag>_50ohm.s2p` | same S-parameters renormalised to a flat 50 Ω |
| `<tag>_sparams.csv` | `freq_hz, S11_mag_db_50ohm, S11_phase_deg_50ohm, S21_..., S11_mag_raw_Zoref` |
| `<tag>_Zo.csv` | `freq_hz, re/im Zpi(port1), re/im Zin=50(1+S11)/(1-S11)` (Zin cols only when the Zo grid matches the S-param grid) |
| `<tag>_convergence.csv` | HFSS adaptive-pass table (pass, tetrahedra, max\|ΔS\|) |
| `<tag>_result.json` | all band-mean numbers + convergence tail |
| `testC_full_design_fit.json` | Test C resonance summary (f_r from eigenmode; why Q_c was not extractable) |
| `testC_eigenmodes.json` | the 5 lowest eigenmodes; mode 1 (8.546 GHz) is the resonator |

## Scripts (run from this directory, in the `qiskit_qc` conda env; they read `../geometry/*.gds`)

| script | what it does |
|---|---|
| `run_ansys_tests.py` | the driver — `python run_ansys_tests.py A B C` |
| `gds_geometry.py` | reads clean metal polygons from the handoff GDS (heals the keyhole ground polygons) |
| `reprocess_zin.py` | recomputes the 50 Ω-referenced `Zin` from a raw `.s2p` + `Zo.csv` without re-solving |
| `testC_eigen.py` | eigenmode solve of Test C -> resonator `f_r` (`testC_eigenmodes.json`) |
| `testC_sweep_only.py` | adds a sweep to Test C's already-solved driven mesh (no re-adapt) |
| `make_summary.py` | regenerates the figures from the result files |

Env knobs: `HFSS_CORES` (default 4), `AEDT_VERSION` (2024.2), `SI_TANDELTA`
(0 for A/B; `1.2e-5` for the Test C resonance run), `HFSS_ZDEF`
(`Zpi`/`Zpv`/`Zvi`), and per-test `TESTA_YHALF`, `TESTC_F0/F1/NPTS/SWEEP`.

## Figures

- `comparison.png` — every method's Z0/Zin on one bar chart
- `testA_line_only.png` — Zpi & Zin vs frequency, plus S-parameters
- `testC_full_design.png` — S21 around the resonance

## Notes / environment quirks

- Driven via `ansys-aedt-core` 1.5 with `settings.use_grpc_api = False` (COM
  transport) — the gRPC transport does not connect on this machine.
- pyaedt's COM release sometimes leaves a stray `ansysedt.exe`; the driver
  kills stray processes and stale `.aedt.lock` files before each run.
- `.aedt` project files: copies in `aedt_projects/`; the originals plus their
  solved `.aedtresults/` folders are in `~/Documents/Ansoft/`.
- Test C's resonator hangs to y ≈ −1317 µm, so its domain uses
  `y_half ≈ 1570 µm` (auto-sized in `run_ansys_tests.py`). `SI_TANDELTA=1.2e-5`
  is set for Test C so the resonator carries a finite Q.
- The Test C driven S-parameters show no resolvable resonance notch; `f_r`
  comes from `testC_eigen.py` (eigenmode). See `SUMMARY.md` for why.
