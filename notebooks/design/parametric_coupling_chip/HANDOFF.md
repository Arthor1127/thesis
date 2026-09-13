# Handoff — Qc extraction for the hanger readout resonator

For the next agent. The previous session ran long; this is the compressed
state plus the single most promising lead, which was found right at the end
and **not yet acted on**: comparing our PALACE config against
**PalaceForCQED** (`/home/ruiz/PalaceForCQED`), a published, hardware-
validated reference workflow the user considers trustworthy.

## Read these first

1. **[`QC_STATUS.md`](QC_STATUS.md)** — the full investigation record (every
   method tried, every bug found, every dead end). Long but complete.
2. **[`README.md`](README.md)** — the directory layout (reorganized
   2026-09-09 into `driven_sweep/`, `eigenmode_sweep/`, `isolation_tests/`,
   `convergence_studies/`, `feedline_sweep/`, `analysis/`, `ansys/`).

## The goal

Map `Qc(coupling_space_um, coupling_length_um)` for the CoupledLineTee-coupled
hanger readout resonator (`readout_resonators.py::build_design`), 42 grid
points, so the chip design can pick a coupling with a target Qc. All
PALACE/FEM solves run on the CNEA SGE cluster via `qsub` — **never locally**
(standing rule, see pinned memory).

## Where it stands

- **`eigenmode_sweep/run_qc_eigenmode.py`** (lossy-lumped-port eigenmode,
  `Qc = fr/(2·Im f)` from the complex eigenfrequency) is the working method.
  Full 42-point map is in **`qc_eigenmode.csv`** (top level), with a
  `ceiling_flag` column.
- **29/42 points (`ok`) are trustworthy** — validated to 1–2% against an
  analytic reference *and* an independent HFSS driven-S11 fit, for Qc up to
  ~2×10⁴, holding to a precisely-bisected ceiling at **Qc ≈ 82,000–85,000**.
- **11 points (`likely_saturated`) + 2 (`borderline`) have NO trustworthy Qc
  from any method we can run.** Above the ceiling the eigenmode `Im f`
  collapses to a numerical noise floor (`loaded_Q` saturates at ~4.5×10⁷
  regardless of the true coupling). The driven-S21 fit never resolved this
  regime either (notches < 0.25 dB, at/below PALACE's absolute S-parameter
  error floor). HFSS, tried both driven and eigenmode by a local ANSYS agent,
  also could not resolve it (`ansys/eigenmode_qc_results/SUMMARY.md`) — its
  adaptive mesh never refines the weak coupler gap.
- The affected design points: `cs22.5_cl150` (976k), `cs25_cl300` (452k),
  `cs15_cl150` (408k), `cs10_cl100/150/250`, `cs25_cl100/150/250`,
  `cs22.5_cl200/250`, and borderline `cs25_cl350` / `cs15_cl250`.

## THE LEAD: match PalaceForCQED's config

`/home/ruiz/PalaceForCQED` — the open-source library behind arXiv:2511.09041
(Ye/Wang/Liu, Tsinghua, "Electromagnetic Feature Extraction in
Superconducting Quantum Circuits ... Using Palace"), benchmarked against real
cryogenic hardware (resonator frequencies within 0.3%, external couplings
within 16% for 3/4 resonators). Julia + Jupyter; runs `palace` directly (no
SQDMetal). Its hand-maintained config templates are in
`examples/defaul_configs/*.json`. **Diff them against what SQDMetal generates
for us** (`out_*/…/outputFiles/config.json` on the cluster, or read
`SQDMetal/SQDMetal/PALACE/Eigenmode_Simulation.py::_create_config_file`).
The differences found so far:

| setting | PalaceForCQED | ours (SQDMetal default / what we ran) |
|---|---|---|
| **`Solver.Order`** | **4** (their convergence sweep only tests orders 3–6) | **1** |
| **`Solver.Linear.Tol`** | **1.0e-8** | 1e-6, relaxed further in some runs |
| **`Solver.Eigenmode.Tol`** | 1.0e-8 | 1e-6 |
| `Solver.Linear.Type` / `KSPType` | `"Default"` / `"GMRES"` | `"SuperLU"` / `"FGMRES"` |
| **`Boundaries.Absorbing`** | **present** (first-order ABC on the air-box outer faces) | **absent** — we use a closed PEC box |
| driven sweep window | ~12 MHz, centred on `fr` (found first by eigenmode), 10 kHz step | 750 MHz, 12001 pts |
| mesh | their validated `r1.5` (`S_min = 1.5·W_tr·2⁻¹·⁵`) | `0.75·½·(smallest feature)` |
| Si `LossTan` | 2.3e-6 (Checchin 2022) | 2.7e-6 |

**Primary hypothesis (testable, cheap):** our Qc ceiling of ~85,000
corresponds to `Im f / Re f ≈ 1/(2·85000) ≈ 6×10⁻⁶` — essentially our
`Linear.Tol = 1e-6`. The solver literally cannot resolve an imaginary
eigenvalue component below its own linear-solve tolerance. **PalaceForCQED
uses 1e-8, two orders tighter, plus order-4 elements.** Either or both should
push the ceiling out enough to reach the 11 stuck points.

**How to test it** — all of these are exposed through SQDMetal, no source
edits needed:
- `PALACE_Eigenmode_Simulation(user_options={"solver_order": 4,
  "solver_tol": 1e-8, "solver_maxits": 300, ...})` — `solver_tol` maps to
  BOTH `Eigenmode.Tol` and `Linear.Tol`; `solver_order` maps to
  `Solver.Order`.
- `eigen.set_farfield('absorbing')` before `prepare_simulation()` — SQDMetal
  writes `Boundaries.Absorbing` (`Model.py:1813`, `:1861`). Currently we use
  the default `'pec'`.
- Order 4 is ~10–50× more expensive (RAM + time) than order 1. Use the
  cluster's big-RAM nodes; start on the cheap isolation test, not the full
  sweep.

**Suggested first job:** rerun `isolation_tests/run_bareres_eigenmode_
isolation.py` (the clean single-port resonator with an analytic
`Qc = (π/4)(R/Z0)`) at `solver_order=4, solver_tol=1e-8`, for a few R values
that gave `Qc ≈ 1×10⁵–1×10⁶` (params in
`isolation_tests/params_bareres_eigenmode_highQ.txt` /
`_bisect2.txt`). If `loaded_Q` now tracks the analytic target past 10⁵
instead of saturating at 4.5×10⁷, the fix is confirmed — then rerun the 11
saturated grid points through `run_qc_eigenmode.py` with the same settings.

Also look at PalaceForCQED's S11→Qc fitting notebooks
(`notebooks/20250406-Palace_simulation_S11_fitting.ipynb`,
`Simulation_Sparameter_fitting.ipynb`, helper `notebooks/mypy.py`) and
`examples/driven.jl` — their driven route (eigenmode to find `fr`, then a
razor-thin driven window + notch fit) is a viable fallback we never tried
properly (we always used a 750 MHz window).

## Confirmed bugs / gotchas from this session (don't rediscover these)

- **SQDMetal units bug**: `retrieve_mode_port_EPR_from_file` returns `kappa`
  in **GHz** despite the docstring saying Hz (Palace's `port-Q.csv` column is
  labelled GHz). `run_qc_eigenmode.py` already multiplies by 1e9 to fix it.
- **SQDMetal `port-Q.csv` needs ≥2 lumped ports** — single-port eigenmode
  jobs get `kappa = NaN` (hardcoded `if len(col_kappa) < 2: raise`), so the
  isolation tests use `eig.csv`'s own `loaded_Q` directly instead.
- **`create_port_2_conds` with `impedance_C` is unreliable** — the config
  JSON gets the right value but the realized capacitance is ~2–3× off even in
  its best-conditioned regime, and unstable/sign-flipping otherwise
  (`isolation_tests/shuntc_validation_design.py`). Don't use it for a
  controlled coupling. The two-port hanger isolation test's Qc numbers are
  invalid because of this.
- **PALACE driven + AMR crashes** (SIGSEGV) in this build (0.16.0), with or
  without our PROM patch — matches open upstream issue awslabs/palace#675.
  PALACE **eigenmode + AMR works** (used in `squid_array_coupling/
  build_and_eigenmode.py`), so coupler-gap-targeted eigenmode AMR is a
  separate untried avenue if the order/tol fix isn't enough.
- **Params-file race condition**: array tasks read their row by
  `SGE_TASK_ID` *when they start*, which can be hours after `qsub`. Never
  edit a params file while its array still has queued tasks — add a new file
  + `PARAMS_FILE=` override instead (every `submit_*.sh` supports this).
- **Never run compute on the cluster login node** — always `qsub`, even for
  a 10-second postprocessing script (`submit_flag_ceiling.sh` etc. exist for
  exactly this).
- The reorganization keeps every `submit_*.sh` doing `cd
  $HOME/circuit_design/parametric_coupling_chip` (top level) and invoking its
  driver by subfolder path (`python -u eigenmode_sweep/run_qc_eigenmode.py`).
  `out_*/` result dirs live at that top level on the cluster and are
  untouched. Local edits → `rsync` the changed subfolder up → `qsub`.

## Refuted hypotheses (don't re-chase)

- "Missing ×2 port-impedance correction" for the CPW lumped port (the
  original `Zin ≈ 40 Ω` vs 50 Ω mystery) — tested `impedance_R` = 25/50/100,
  R=50 was best. HFSS independently gets ~50 Ω, so PALACE's driven lumped-
  port `Zin` is just an unreliable diagnostic here; `fr`/`Qc` don't depend on
  it.
- "Two simultaneous lossy ports break the eigenmode extraction" — tested by
  forcing one port near-open; the `loaded_Q` vs `port-Q.csv` discrepancy
  persisted. The mode-hybridization explanation (a competing eigenmode
  ~150–250 MHz away) has moderate supporting evidence
  (`analysis/mode_gap_vs_ratio.csv`, corr ≈ −0.49) but is unconfirmed.
- "A simple analytic `Qc(coupling_length)` model to fit-and-extrapolate" —
  tried a weak-coupling power law and a junction+distributed model
  (`analysis/plot_qc_length_fits.py`, `images/qc_length_fits.png`). Neither
  fits (R² ≈ 0.15). Our grid confounds `coupling_length` with the coupler's
  `down_length` (held to a fixed 700 µm budget, `corr = −1.000`), so we've
  never varied coupling length in isolation.

## ANSYS side (external agent, user relays)

- `../ansys_impedance_handoff/EIGENMODE_QC_REQUEST.md` Part 2 (lossy-port
  eigenmode on `cs22p5_cl250` / `cs25_cl100`) — **answered, negative**, see
  `ansys/eigenmode_qc_results/`. HFSS can't resolve it either on their 16 GB
  machine.
- Nothing further is pending on the ANSYS side unless you send a new request
  (they have `run_ansys_tests.py` + the GDS geometries and can rerun with new
  settings — e.g. their own `RESUME.md` proposes a coupler-seeded eigenmode
  on ≥24 GB, or a driven-S21 + coupler-seed + narrow sweep).
