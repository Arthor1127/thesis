# Status — schizo-resonator grounding-position sweep

Last updated: end of the session covering install through mesh-memory
debugging. Update this file at the end of each future session so the next
one starts from an accurate picture rather than re-deriving it.

## Goal

For a single-strap, movable-ground-position half-wave resonator (both
ends open, `otg1`/`otg2`), produce:
1. The baseline (ungrounded) S21 spectrum.
2. A colormap: x = grounding strap position along the resonator,
   y = frequency, z = |S21| — showing how the spectrum evolves as the
   strap moves, across a fixed 4.5–20.5 GHz window covering the 1st–4th
   harmonics simultaneously for every position.

Pipeline: `build_design.py` → `build_and_eigenmode.py` (Stage A, per-
position eigenmode reference) → `build_and_driven.py` (Stage B, per-
position driven S21 sweep) → `aggregate_and_plot.py` (Stage C, final
plots). SGE array wrappers: `stageA_eigenmode.sh`, `stageB_driven.sh`,
`stageC_aggregate.sh`, driven by `run_full_sweep.sh`.

## Calibration — DONE, confirmed good

`TOTAL_LENGTH_MM = 14.9669` in `build_design.py` is confirmed correct:
baseline (`ground_pos=None`) eigenmode run measured `f1 = 5.0949 GHz`
vs. 5.0000 GHz target (**+1.90%**, well inside the 5% acceptance
threshold). Do not re-derive this unless the design itself changes
(trace width/gap, substrate, chip size, etc.).

## Mesh settings — validated on the eigenmode side, driven side in progress

Current defaults in both `build_and_eigenmode.py` and
`build_and_driven.py`:
- `--bg-min-size-um 15 --bg-max-size-um 150` — background mesh across
  the whole domain. Validated via isolated `gmsh`-only testing (see
  LESSONS_LEARNED.md): 84% fewer elements than the original 8/100um
  setting, 6.3x faster mesh generation, with `fillet_resolution` making
  no measurable difference either way.
- `--strap-margin-um 3.0` — **fixed** margin around the strap's
  bounding box for its local fine-mesh region. This replaced an earlier
  clearance-derived version that had a critical bug (collapsed to 0
  margin when the measured "clearance" was to the strap's own
  attachment point, causing an instant mesh-validity crash). Do not
  reintroduce clearance-based auto-margin logic without addressing that
  flaw properly (would need to exclude the local trace segment under
  the strap itself from any distance calculation, not just take
  nearest-distance-to-own-component blindly).
- `--resonator-min-size-um 2.0` (driven) / hardcoded `2e-6` (eigenmode)
  — fine mesh size inside the strap-local box itself.
- `--solver-order 2` — kept at user's request (declined to downgrade
  for memory reasons). This is the most memory-expensive setting still
  in play; see below.

## Driven simulation memory — NOT YET CONFIRMED WORKING END TO END

As of the last real (pasted, verified) log:
- The corrected mesh settings above successfully produced a well-formed
  mesh (1660 boundary elements, matching the one known-good reference
  run) rather than the earlier 42K-97K-element blowups.
- `Assembling system matrices` numbers with the corrected mesh have
  **not yet been directly observed and confirmed** in a pasted log — a
  prior message in this session incorrectly stated specific numbers
  (H1=1.47M etc.) without them having actually been pasted; those
  figures were fabricated/predicted, not observed, and should be
  disregarded. Do not trust any memory/unknown-count figures for the
  driven sim until they appear in an actual pasted log.
- No successful driven-sim run (i.e., one that reaches
  `retrieve_data()` and produces real S11/S21 numbers) has occurred yet
  in this session, at any mesh/node combination tried so far.
- Test job in flight at time of writing: `test_driven8`
  (`test_single_driven.sh`, `--ground-pos 5.906 --solver-order 2`,
  no node pin, `~/sweep_v2_test/driven_5.906_v8`). Check its outcome
  first before drawing any conclusions or proceeding to the full array.

## Node notes

- Avoid `compute-4-7`, `compute-4-8` (small RAM, ~31G) for anything
  beyond trivial jobs.
- `compute-4-9` had an `E` (error) state flag as of last check — avoid
  until confirmed healthy again.
- `compute-4-18` (503G) is the largest node and was confirmed idle
  multiple times, but a job pinned specifically to it once stalled in
  `qw` for 13+ hours despite the node showing genuinely idle the whole
  time (`qhost`/`qstat -f` both clean) — a likely scheduler-side issue,
  not a resource conflict. If this recurs, don't keep re-pinning to a
  single host; let SGE choose freely, or escalate to
  `soportefisica@cnea.gob.ar` with the specific evidence (job ID,
  wait time, confirmed node idleness).

## Sweep parameters (once the driven sim is confirmed working)

- `sweep_positions_mm(n_points=9)` — evenly spaced across
  `[0.5mm, TOTAL_LENGTH_MM - 0.5mm]`, i.e. currently
  `[0.5mm, 14.4669mm]`. Plus the `none` (ungrounded) baseline = 10
  positions total for Stage A.
- Stage B (driven): one continuous 4.5–20.5 GHz sweep per position
  (not split by harmonic anymore — an earlier per-harmonic-chunk design
  was replaced), currently 321 points (~50 MHz spacing) — chosen as a
  starting point, **not verified against the driven sim's actual loaded-
  Q linewidth** (the eigenmode sim's unloaded Q ≈ 4×10^5 would need
  ~12.5 kHz spacing to fully resolve, computationally infeasible at this
  span; the real driven/loaded linewidth is expected to be broader but
  is unconfirmed). If a completed driven sweep shows no visible peaks,
  increase `--n-points` or narrow the range around empirically-found
  peak locations.
- 10 positions × 1 continuous sweep = 10 Stage B tasks (not 40 — an
  earlier 4-harmonics-as-separate-chunks design was replaced with one
  continuous sweep per position).

## Immediate next steps, in order

1. Check the outcome of `test_driven8` (or whatever the latest single-
   position test job is) — get the actual pasted log, not an assumption.
2. If it succeeds: inspect the real S21 result for sanity (does it show
   plausible resonance features anywhere in 4.5-20.5 GHz?), then note
   the *actual* peak memory reported (`Estimated peak per-rank/per-node
   memory usage`) and use that to decide which nodes the full 10-task
   Stage B array can safely target (aim for headroom on more than just
   the single largest node, so tasks can run in parallel rather than
   serializing on one host).
3. If it fails again: read the actual failure mode carefully (OOM vs.
   mesh-validity vs. something new) before changing anything — don't
   assume it's the same bug as last time.
4. Once one driven position is confirmed working end to end: submit
   the real `stageA_eigenmode.sh` + `stageB_driven.sh` (10 tasks each)
   with the validated node list, then `stageC_aggregate.sh` for the
   final plots.
5. Inspect `~/sweep_v2/plots/baseline_S21.png` and
   `S21_evolution_colormap.png` for the actual physics result.
