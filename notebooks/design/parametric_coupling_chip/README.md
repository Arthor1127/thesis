# parametric_coupling_chip

Hanger-resonator Qc(coupling_space, coupling_length) extraction for the
readout-resonator design — PALACE/SQDMetal simulations on the CNEA cluster,
cross-checked against ANSYS/HFSS. **Start with [`QC_STATUS.md`](QC_STATUS.md)**
for the full narrative (what was tried, what worked, what's still open).

## Layout

Every experiment family lives in its own subfolder: `run_*.py` (the PALACE
driver, submitted as a cluster job), `submit_*.sh` (its `qsub` script),
`params_*.txt` (point lists for array jobs), and any `aggregate_*.py` /
`refit_*.py` / `flag_*.py` postprocessing that turns raw per-point results
into a summary CSV or figure.

| Folder | What it is |
|---|---|
| `driven_sweep/` | The original driven-S21 notch-fit method (`run_qc_sim.py`) and its production sweep — works, but hits a hard shallow-notch reliability ceiling on this design (see `QC_STATUS.md`). Includes `check_line_impedance.py` (bare-CPW-line Zin isolation test, unrelated to Qc). |
| `eigenmode_sweep/` | The lossy-eigenmode method (`run_qc_eigenmode.py`) — the one we trust, and the production 42-point Qc map (`qc_eigenmode.csv`, top-level). |
| `isolation_tests/` | Clean analytic-ground-truth validation geometries: a bare single-port resonator (`run_bareres_eigenmode_isolation.py`), a two-port hanger with a lumped capacitor (`hanger_isolation_design.py` — numbers now known invalid, see `QC_STATUS.md`), and a resonance-free shunt-capacitor check (`shuntc_validation_design.py`) that confirmed why. |
| `convergence_studies/` | Solver-setting and mesh-resolution convergence sweeps that produced the validated production recipe (`mesh_factor=0.75`, relaxed tolerances). Reuses `driven_sweep/run_qc_sim.py` and `driven_sweep/check_line_impedance.py` as its underlying drivers. |
| `feedline_sweep/` | An earlier eigenmode cross-check (no ports, closed-box) used to pin the bare resonator frequency; superseded by the isolation tests above. Includes raw archived PALACE output (`output/`, `out_fdln/`). |
| `analysis/` | Post-hoc diagnostics: mode-gap-vs-driven/eigenmode-disagreement correlation, and the analytic-model fit attempts (both found not to fit — see `QC_STATUS.md`). |
| `ansys/` | Everything delivered by the local ANSYS/HFSS cross-check agent: the original impedance handoff results, the bare-resonator S11 cross-check, and the (negative but conclusive) lossy-eigenmode cross-check. |
| `images/` | Every figure referenced from `QC_STATUS.md` — all plotting scripts save here directly. |

Top-level: `QC_STATUS.md`/`.pdf` (the report), `readout_resonators.py`/`.ipynb`
(the shared hanger-resonator geometry builder, imported by scripts across
several subfolders — kept at the root for that reason), `full_chip.ipynb`,
`qc_eigenmode.csv` (the current best Qc(coupling_space, coupling_length) data
product, cross-cutting enough to stay at the root), and `images/`.

The `notebooks/design/ansys_impedance_handoff/` and
`notebooks/design/ansys_eigenmode_s21_crosscheck/` handoff directories (with
the outgoing `HANDOFF.md`/`EIGENMODE_QC_REQUEST.md` requests and exported GDS
geometry for the ANSYS-side agent) are siblings of this whole directory, not
inside it.

## Running things

Every `submit_*.sh` still assumes it's `qsub`'d from the cluster's
`~/circuit_design/parametric_coupling_chip/` (the top level of this same
tree, mirrored via `rsync`) — they `cd` there internally regardless of which
subfolder the script itself lives in, then invoke their Python driver and
params file by that driver's actual subfolder path (e.g.
`python -u eigenmode_sweep/run_qc_eigenmode.py`). Existing `out_*/` result
directories (e.g. `out_qc_eig/`, `out_qc_prod/`) live at that same top level
on the cluster and are untouched by this reorganization — they're not synced
back locally (too large), only their aggregated CSV/JSON/figure outputs are.

Standing rule: PALACE/FEM/MPI solves always run on the CNEA SGE cluster via
`qsub`, never locally — see the pinned memory on this. Geometry/mesh
inspection and CSV postprocessing are fine locally.
