# `Analisis_gaps_cris` stanene gap map — SGE cluster pipeline

Runs the (B_x, ε_site) indirect-gap map from `../Analisis_gaps_cris.ipynb`
(stanene s+p Slater–Koster model) on the CNEA SGE cluster.

## Layout

| file | role |
|------|------|
| `gaps_kernel.py`  | `stanene_pythtb`, `gap_2d`, `gap_1d` — copied verbatim from the notebook |
| `build_grid.py`   | writes `grid.json` (mesh 50, B∈[0,1], ε∈[0,0.6], mallado 240, ribbon 35 cells) |
| `gaps_sweep.py`   | one array task = one B row (50 ε points), checkpointed → `results_stanene/row_III.npz` |
| `make_plots.py`   | assembles rows → `gap_map{1D,2D}.npy`, `B_vals.npy`, `E_vals.npy`, `gap_map.npz`, thesis-styled `gap_map_{1D,2D}.{png,pdf}` |
| `submit_gaps.sh`  | qsub the sweep array (`-tc 96`, 1 slot/task, `mem=4G`) + a `-hold_jid` plot job |
| `deploy.sh`       | rsync this dir to `cnea-cluster:~/unitary/cris/` and optionally submit |

## Run

```bash
./deploy.sh submit
```

Watch: `ssh cnea-cluster 'qstat -u $USER'`; done when
`ls results_stanene/row_*.npz | wc -l == 50`. Results + figures land in
`results_stanene/` on the cluster.

## Notes

- The `cris-gaps` conda env (pythtb 2.0 + numpy/scipy/joblib/matplotlib +
  `thesis_style`) already exists on the SGE cluster — no setup job needed.
- 1 slot per task: small dense linalg, single-threaded BLAS is fastest and
  1-slot tasks backfill with near-zero queue wait. `highmem` excluded (other
  agent's).
- `gap_2d` is the **indirect** gap `min(cond) − max(val)` and is stored raw —
  it can be negative (band overlap / metallic); the plot clips at `vmin=0`.
- Per-task peak RAM ≈ 3–3.5 GB (the `(240, 560, 560)` ribbon Hamiltonian batch).
- The earlier germanene run's results are in `results/` on the cluster
  (`Analisis_gaps.ipynb`, mesh 100) — kept separate.
