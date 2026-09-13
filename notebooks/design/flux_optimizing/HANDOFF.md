# Handoff — SQUID-array inductive coupler, flux-threading sweep + corroboration

Status as of **2026-09-04 21:55 -03**. Everything below is reproducible from
this directory (`notebooks/design/flux_optimizing/`) plus the mirrored copy on
the cluster (`arturo.ruiz@10.73.25.223:circuit_design/flux_optimizing/`).

**Update:** this project is being closed out at the state described below.
The magnetostatic sweep is complete for practical purposes (2400/2400 grid
points, 0 failures, 46 flagged not-yet-reconverged — a small fraction whose
correction is expected to be a near-uniform multiplier per §4, so it doesn't
change the ranking). The remaining job chain (§2, as it stood) was cancelled
by user decision rather than run to completion: the 3 still-running gap-fill
tasks were projected at 3–7h to finish (measured GMRES rate ~0.05–0.1
iter/s, hundreds of iterations short of the 2000-iteration cap) and the
downstream reconverge + real-design-eigenmode corroboration stages were
judged no longer worth the cluster resources. Local `sweep_results.csv` /
`sweep_colormaps.png` were re-verified against a fresh pull of all 300 shard
CSVs and are byte-for-byte consistent with this final state — no re-fetch
needed.

## 1. What this project is

Optimizing a 23-loop `SquidArray` inductive coupler read out by a hairpin CPW
resonator. Two parallel efforts:

1. **5-D magnetostatic sweep** (`sweep_flux.py`) — a uniform-current
   magnetostatic Palace solve over `route_clearance, segment_length,
   squid_gap, arm_semi_length, arm_width` (4×5×6×5×4 = 2400 points), ranking
   geometries by `M = ΣΦᵢ/I` (`Phi0_per_mA` column). Array metal is PEC, loops
   left **open** (no Josephson junctions — negligible screening at the design
   `L_J`, see §4).
2. **Eigenmode corroboration** (`run_eigenmode.py` / `extract_eigen_M.py` /
   `run_eigenmode_realdesign.py`) — same open-loop geometry, but the resonator's
   own self-consistent standing-wave current replaces the uniform-current
   assumption. Extracts `M = ΣΦᵢ/I_short` with the reference current sampled at
   the **shorted end** (user's explicit choice — comparable across loops, unlike
   a per-loop local reference).

## 2. Cluster job chain — CANCELLED, closed out

```
1938674 (main sweep, 300 tasks)          -- DONE
   -> 1938675 (gap-fill pass 1, 150 tasks) -- CANCELLED 21:49 (3 tasks still running:
                                              004, 149, 150; measured ~0.05-0.1 iter/s,
                                              projected 3-7h to finish — not worth it)
      -> 1938676 (gap-fill pass 2, 150 tasks, held)     -- CANCELLED, never started
         -> 1938677 (aggregate, held)                    -- CANCELLED, never started
      -> 1938713 (reconverge 46 stragglers, 25 tasks, held) -- CANCELLED, never started
            -> 1938714 (final aggregate, held)           -- CANCELLED, never started
               -> 1938727 (real-design eigenmode, held on 1938714) -- CANCELLED, never started
```

Also cancelled: `1938728` (task 141 requeue after a `be_copahue` suspension).

No further cluster work is planned for this project. If it's ever revisited,
the 46 not-yet-reconverged points and the real-design eigenmode corroboration
(`run_eigenmode_realdesign.py`, still fully written and usable) are where it
left off — see §5.

## 3. Files (all in this directory unless noted)

Build / geometry:
- `build_design.py` — the design builder. `geometry_anchors()` (array/hairpin
  layout, added this session) and `loop_aperture_polys()` are the two reusable
  geometry helpers. `array_ground_margin_um` (new default 6 µm rectangular
  ground pocket) fixes a real bug (SquidArray's own conformal pocket pinches
  in between loops and floods grounded tabs there) — **the main sweep is
  deliberately pinned to the OLD conformal pocket** (`sweep_flux.py` passes
  `array_ground_margin_um=None`) for internal consistency across all 2400
  points; only new/future runs get the fixed rectangular pocket.

Magnetostatic sweep:
- `sweep_flux.py`, `submit_sweep.sh` / `submit_sweep_hi.sh` (gap-fill) —
  the sweep itself.
- `reconverge_stragglers.py`, `submit_reconverge.sh` — re-solves the 46
  points that hit the old `SOLVER_MAXITS=1200` cap (before it was bumped to
  2000 mid-sweep) with a boosted `solver_maxits=6000, gmres_maxsize=300`.
- `aggregate_sweep.py`, `submit_aggregate.sh` — merges shards →
  `sweep/sweep_results.csv` + `sweep/sweep_colormaps.png`. **Now uses
  `thesis_cmap`** (fixed this session).
- `viz_sweep_final.py` — local trends/fits figure from `sweep_results.csv`,
  thesis-styled. Run as `python viz_sweep_final.py sweep_results.csv`.
- `count_done.py` — quick done/error tally across shards.

Geometry gallery (all 2400 geometries rendered + sanity-checked, **0 flagged**):
- `render_geometry.py`, `submit_render.sh`, `submit_sheets.sh` —
  `geometry_gallery/_sheets/clr*_aw*.png` (16 contact sheets). **Not yet in
  thesis style** (categorical material-identity diagram, not a data plot —
  flagged as skipped, not forgotten).

Eigenmode corroboration:
- `run_eigenmode.py` + `extract_eigen_M.py` — one flux_optimizing sweep-corner
  eigenmode run + B-field extraction. Both thesis-styled now.
- `run_eigenmode_realdesign.py` — same method, but on `squid_array_coupling`'s
  actual target design (N=23, squid_gap=25um, meandered readout,
  `RES_LENGTH_MM=3.74575` → their own converged `f_r=7.9355 GHz`). Thesis-styled.
  **This is `1938727`, still queued.**
- `submit_eigenmode.sh`, `submit_eigenmode_realdesign.sh`, `submit_extract.sh` —
  submit scripts. Current default: **32 cores / 160 G** (bump to 220G, which
  forces the single `highmem` node, for configs at/near the biggest-loop
  corner — see the comment block in `submit_eigenmode.sh`, sized from a
  measured 183 G peak at the `best` config).
- `replot_eigen_M.py` — redraws a saved `eigen_M_cyc*.npz` locally in thesis
  style, no cluster round-trip. Used for `eigen/small` and `eigen/best`.

Results already fetched locally (safe to work with now):
- `sweep_results.csv`, `sweep_colormaps.png` — **snapshot from the current
  2400/2400-attempted, 46-not-yet-reconverged state** (re-fetch once
  `sweep/DONE.txt` appears on the cluster for the final version).
- `sweep_trends_final.png` — trends/fits figure from that snapshot.
- `eigen/small/`, `eigen/best/` — completed eigenmode corroboration runs (see §5).
- `eigen/realdesign/` — **not yet populated**, waiting on `1938727`.

## 4. Key results so far

**Scaling law** (full 2400/2400 attempted, log-space R²=0.99, 2% median err):

```
Phi0/mA ~ 6.74e3 * route_clearance^-2.95 * squid_gap^+0.64 * arm_semi^+0.72
                 * arm_width^+0.27 * (seg+3)^-0.04
```

`route_clearance` (wire↔array gap, ∝1/r³ dipole falloff) dominates by far;
`squid_gap`/`arm_semi` (loop area) show diminishing returns (~0.65–0.7, not
linear); `arm_width` (flux-funnelling) is weak and saturating (~0.27);
`segment_length` is essentially flat/negligible.

**Best configs** (`clr30/.../gap40/semi24/aw8` corner):
- Best **converged**: `seg=3` → 47.49 Φ0/mA = **98.2 pH**
- Best **overall** (not yet reconverged, pending `1938713`): `seg=0` → 49.53
  Φ0/mA = 102.4 pH

**Eigenmode corroboration (current-profile-weighting correction):**

| config | M (magnetostatic) | M (eigenmode, shorted-end ref) | ratio |
|---|---|---|---|
| smallest loops | 7.28 pH | 4.36 pH | 0.598 |
| best-M config | 98.20 pH | 57.13 pH | 0.582 |

Two very different geometries land within 3% of the same ~0.58–0.60
correction factor — good sign it's dominated by the shared elbow/short layout,
not loop-geometry specifics. **Practical read: best-guess real M ≈ 57 pH**, not
98 pH — but the *ranking* from the magnetostatic sweep should carry over
essentially unchanged (a near-uniform multiplier doesn't reorder configs).
Per-loop flux is NOT a simple monotonic short→open decay — it's highest near
the grounded elbow end, dips to a broad minimum ~loop 19–21, ticks back up near
the terminations (folded-resonator current profile, not a straight quarter-wave
line). Still waiting on: the `squid_array_coupling` real-design corroboration
(`1938727`) for a third, independent data point at the design's actual
operating parameters.

**Auto-pick bug found and fixed:** for the `best` config, picking "the mode
with the most array-flux-linkage" chose a spurious *array self-resonance*
(29.35 GHz — the array's own leads+ground-capacitance ladder, matching a
`squid_array_coupling` result: `L_J=15nH → phonon band top ~30 GHz`) instead of
the true resonator fundamental (8.74 GHz, matching the analytic bare-λ/4
estimate almost exactly). `extract_eigen_M.py` and `run_eigenmode_realdesign.py`
now pick by frequency-proximity to the analytic estimate instead, with a
printed audit trail of both criteria.

## 5. Open questions / natural next steps

1. ~~Finish the sweep chain~~ — **closed out, not pursued.** Cancelled by
   user decision on 2026-09-04 (§2): the 46 not-yet-reconverged points and the
   remaining gap-fill/reconverge stages were judged not worth further cluster
   time given the ranking is already stable. `sweep_results.csv` /
   `sweep_colormaps.png` reflect the final state (2400/2400, 0 failures, 46
   flagged) and need no further fetch.
2. ~~`1938727` real-design eigenmode~~ — **closed out, not pursued.**
   Cancelled before it ever ran. If picked back up later,
   `run_eigenmode_realdesign.py` is still written and ready — their array
   params (squid_gap=25, arm_width=4, arm_semi=15, segment_length=6) are the
   flux_optimizing module *defaults*, so `loop_aperture_polys()` with those
   values reproduces the right apertures for any manual cross-check.
3. **Refresh `MODEL.md`** — several numbers in it are still marked
   `[partial]`; can now be finalized against the closed-out sweep state above
   rather than waiting on a `DONE.txt` that will never appear.
5. **Kinetic inductance follow-up (discussed, not started):** `Model.py` has
   `add_kinetic_inductance(Ls)` (Palace's native `Impedance` boundary, `Ls` in
   H/□ — equivalent to specifying the Pearl length, `Λ=2·Ls/μ₀`, up to the
   thin-film/local-boundary caveat discussed in chat). `SQDMetal.Utilities.
   Kinetic_Inductance` has a full Mattis-Bardeen calculator to derive `Ls` from
   material parameters. Not yet wired into any script here — would need a
   real `Ls` value for the actual film (grAl — see the granular-aluminum
   memory) and a paired PEC-vs-impedance run to extract the kinetic-inductance
   fraction `α = (L_total-L_geo)/L_total`.
6. **Φ₀/2 flux-bias / array self-resonance risk** (discussed, not
   implemented): biasing near Φ₀/2 raises the effective Josephson inductance a
   lot, which (per the `L_J=15nH→30GHz, 600nH→4.7GHz` scaling already
   established in `squid_array_coupling`) will sweep the array's own
   self-resonance down through the readout frequency. Real design risk worth
   flagging when real junction parameters are chosen; not something to act on
   yet.
7. `render_geometry.py` gallery diagrams could be converted to the thesis
   categorical palette (`ts.CYCLE`) if wanted — flagged, not done.

## 6. Standing rules (don't relitigate these)

- **Never run simulations locally** — every Palace/FEM solve goes through
  `qsub` on the cluster, even small validation ones. Geometry build, CSV
  post-processing, and plotting are fine locally.
- **Never work directly on the cluster frontend** — always via a submitted
  job, even for package installs.
- **All figures use `thesis_style`/`thesis_cmaps`** (`directives/figure-style.md`)
  — `ts.use()`, `cmap='thesis_cmap'`, no hand-picked colors, no `viridis`.
  Confirmed installed in the cluster's `qcg-quantum-design` conda env too.
- Default cluster sizing: 32 cores for a single solve (this account's standing
  default), small-core jobs fanned across the ~200-slot quota for many-point
  sweeps. Size memory from measured evidence, not guesswork.
