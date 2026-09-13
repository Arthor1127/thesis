# Handoff: reproduce PalaceForCQED's driven S11 plot in SQDMetal

## Goal

`/home/ruiz/PalaceForCQED` (hardware-validated, arXiv:2511.09041) has a real,
published driven-S11 result — `examples/imgs/driven_new.svg`, several clean
resonance dips at -6 to -9dB in `|S[1][1]|`. Build the **same design** in
SQDMetal/qiskit_metal and reproduce that same S11 curve. If SQDMetal can match
a design PalaceForCQED has already validated against real hardware, that
becomes a second, independent, real-device-anchored cross-check — on top of
the from-scratch Julia/raw-Gmsh line crosscheck this session already
completed (see `QC_STATUS.md`'s "Update (2026-09-12)" section). The point is
to **rule out knobs**, not add more: a design with a known-good published
reference removes the "which parameter is actually right?" ambiguity a
synthetic test always carries.

## The PalaceForCQED reference is fully available — geometry included

Confirmed by reading the source directly (not guessing). Original location:
`/home/ruiz/PalaceForCQED` (a separate repo on this machine). **Local copies
of every file below are already sitting in `palace_cqed_reference/` in this
same directory**, so this task is self-contained and doesn't require that
other repo to be present/synced — paths given here are the local copies;
re-fetch from the original repo only if you need something not copied (e.g.
`Impedance_based_analysis.ipynb`, or to confirm something against the live
git history).

- **Geometry**: `palace_cqed_reference/2DQv8_eb4_edit.geo` (originally
  `examples/geo_files/2DQv8_eb4_edit.geo`) — a real CAD/GDS-style geometry
  file (20,528 lines), loaded by `MeshGenerator.jl`'s `generate_mesh`
  (`palace_cqed_reference/MeshGenerator.jl`, originally `src/MeshGenerator.jl`
  — copied since `generate_mesh` is the function that actually consumes the
  `.geo` file). The device is defined parametrically in
  `palace_cqed_reference/def_2DQv8_eb4.jl` (originally
  `examples/def_2DQv8_eb4.jl`): it's a full transmon-qubit + readout-resonator
  + coupling structure, with named regions (`RQ_RECT`, `READOUT_RESONATOR_RECT`,
  `RESONATOR_RECT`, `QUBIT_RECT`, `JUNCTIONS`, `EXCITAION_LUMPED_PORT`, a wave
  port at the domain boundary `WAVE_RECT`) — not a bare CPW line like this
  session's crosscheck.
- **Solver config**: `palace_cqed_reference/2DQv9_eb4_driven_lumped.json`
  (originally `examples/defaul_configs/2DQv9_eb4_driven_lumped.json`; note:
  v9-named config file reused against the v8 geometry — that's how their own
  `run.jl` calls it, not a typo to fix).
- **Run recipe**: `palace_cqed_reference/run.jl` (originally
  `examples/run.jl`), function `driven_lumped(...)` — this is the literal
  function that builds the mesh (`generate_mesh(mesh_config)` with
  `GEO_PATH` from `def_2DQv8_eb4.jl`) and drives the sweep. Read this
  function fully before starting; it has the exact frequency-window/adaptive
  settings used for the published plot (several narrow windows are
  commented out at the bottom of the file, e.g. `driven_lumped(7.0, 7.8,
  0.001, 0, 16.78543, 21.70551)` — these look like the actual per-dip
  narrow-window recipe, worth checking against the multiple dips visible in
  `driven_new.svg`).
- **Plot**: `palace_cqed_reference/driven_plot.jl` (originally
  `examples/driven_plot.jl`) — confirmed by the exact title-string match
  (`"Port Response of 2DQv8_eb4"`) that this is what produced
  `driven_new.svg`. Reads a `port-S.csv`-style CSV via CSV.jl/DataFrames.
- **The reference image itself**: `palace_cqed_reference/driven_new.svg`
  (originally `examples/imgs/driven_new.svg`) — the actual target to match.

**Key parameter already extracted**: `DEFAULT_ORDER = 4` in
`def_2DQv8_eb4.jl`. PalaceForCQED runs this reference design at **solver
order 4**, not 1 or 2. That's directly relevant — see "What we already know"
below.

Not yet checked: the exact `Absorbing`/`PEC` boundary setup for this specific
design (the earlier session found PalaceForCQED's *eigenmode* configs
explicitly comment out absorbing boundaries but its *driven* configs use
them — verify which applies here by reading the actual generated config, not
assuming), and whether `MeshGenerator.jl`'s mesh-sizing approach differs from
`impedance_crosscheck/mesh_si_cpw_line.jl`'s Distance+Threshold field pattern
in the same way SQDMetal's did (see below).

## What this session already established (read `QC_STATUS.md`'s "Update
(2026-09-12)" section in full before starting — this summary is not a
substitute)

A parallel investigation this session built a much simpler bare-CPW
through-line in both a hand-written raw-Gmsh/Julia model
(`impedance_crosscheck/`) and SQDMetal (`driven_sweep/check_line_impedance.py`,
`MATCH_JULIA=1` mode) and got them to agree. Findings directly relevant here:

1. **Solver order is the single biggest lever for SQDMetal accuracy on
   impedance/S-parameter quantities.** Order 1→2 alone closed ~15 percentage
   points of a ~20% Z0 error, on an otherwise-identical, fully mesh-converged
   model. PalaceForCQED's own reference design runs at **order 4** — expect to
   need at least order 2, quite possibly order 4, to reproduce their S11
   curve. Order 4 costs roughly 15x order 2 per solve on the same mesh (by
   extrapolation from the confirmed 1→2 cost ratio of ~14.9x) — **budget
   cluster time accordingly, and size the job's memory request generously**
   (see `cluster-resource-request-sizing` project memory).
2. **`set_farfield("absorbing")` is never called by default anywhere in this
   codebase's SQDMetal usage** (`Frequency_Driven_Simulation.__init__` sets
   `self._ff_type = {}`) — every outer boundary defaults to PEC. This made no
   difference in the simple bare-line test (walls were already far enough
   away), but a real qubit+resonator+readout geometry has near-fields at
   different scales — **do not assume this is harmless here without
   checking**, the way it was safely ruled out for the bare line.
3. **`Z_PROP_BOTTOM` defaults to `0.0`** in this project's SQDMetal usage
   pattern (`set_zBoundary_as_proportion(top, bottom)`), meaning the domain
   floor sits directly at the substrate's back face with no air gap unless
   set otherwise. Check whether the real device is meant to have a grounded
   backplane or an open substrate before deciding whether this matters.
4. **Port fin length** (`create_port_CPW_on_Route`/`_on_Launcher`'s
   `len_launch`) should be kept thin relative to the actual CPW gap it's
   filling, per that function's own docstring — the old default (30µm) was 5x
   a 6µm gap in the bare-line test and was backwards from the guidance, though
   its effect there was small.
5. **Mesh density/taper parameters were NOT the source of the bare-line's
   discrepancy** — near-line size, taper width, and far-field size were all
   tuned independently; none moved the impedance result by more than ~2%.
   Don't spend much time here unless order/boundary/ground-plane-width are
   already ruled out. **Separately — adopt as a standing practice regardless
   of whether it matters for this task's outcome**: Julia's `Distance` field
   (`mesh_si_cpw_line.jl`) sources from the gap surfaces specifically
   (narrow, localized) and combines `PointsList`+`CurvesList`+`SurfacesList`
   jointly, with `Sampling=ceil(L/l_trace)` scaled to the target resolution
   (~530 samples for a 4mm line). SQDMetal's `fine_mesh_components` sources
   from the whole named component's footprint (wider) and only ever sets one
   of `SurfacesList`/`CurvesList` (never combined with points), with a flat
   `Sampling` default of 150 regardless of geometry size
   (`gmsh_dist_func_discretisation`, now env-overridable as
   `GMSH_DIST_SAMPLING` in `check_line_impedance.py`) — confirmed by direct
   node-density cross-section comparison to produce a visibly
   blockier/rectangular fine-mesh region than Julia's smooth radial cone.
   For the real qubit+resonator geometry in this task, prefer sourcing fine
   meshing from the narrowest relevant feature (the actual coupling
   gap/junction regions, not whole components) and scale `Sampling`
   accordingly, even if — as in the bare-line test — it doesn't turn out to
   be what explains any mismatch you find.
6. **When comparing S-parameters to a Z0/impedance target, use a proper
   two-port NRW line inversion** (`impedance_crosscheck/stage2_analyze.py::
   invert_line`) if the reference design allows it (a uniform through-line
   does; a resonator+qubit design generally does not — for this task you're
   comparing `|S11|(f)` shapes/dip depths/frequencies directly against
   `driven_new.svg`, not inverting to a single Z0 number).
7. **Both production Qc scripts** (`driven_sweep/run_qc_sim.py`,
   `eigenmode_sweep/run_qc_eigenmode.py`) default to `solver_order=1` and
   never call `set_farfield` — flagged but **not yet fixed** in this
   session. Out of scope for this handoff's task but worth knowing since
   whatever you learn here about order/boundary sensitivity for a
   qubit+resonator geometry (closer to the real production design than the
   bare line was) is directly relevant feedback for that separate,
   still-open fix.
8. **No SQDMetal source was modified this session.** Every fix was a
   usage-pattern change (which `user_options`/API calls get made). Keep it
   that way unless you find a genuine SQDMetal bug with equally direct
   evidence (this session found port coordinate math, fin area, materials,
   and LumpedPort JSON structure all correct by direct empirical
   verification against generated files — check real output before
   suspecting the library).

## Suggested approach

1. Read `examples/run.jl`'s `driven_lumped` function and `def_2DQv8_eb4.jl`
   fully. Understand the real geometry (transmon capacitor pads, junction
   rectangles, readout resonator, coupling structure, lumped excitation
   port) well enough to redraw it as qiskit_metal components — this is the
   main work of the task, likely a `TransmonPocket`-family qubit +
   `RouteMeander`/`RoutePathfinder` resonator + a coupling structure, sized
   to match the `.geo` file's actual dimensions (extract them directly from
   the `.geo` file or from the `Rectangle(...)` calls in
   `def_2DQv8_eb4.jl` — those are already real coordinates in µm).
2. Build the equivalent SQDMetal `PALACE_Driven_Simulation`, starting from
   whichever of `run_qc_sim.py` / `check_line_impedance.py` is the closer
   template (probably `run_qc_sim.py`, since it already builds a real
   resonator+coupler, not just a bare line).
3. Apply the guidelines above proactively rather than rediscovering them:
   start at `solver_order=2` minimum (budget for testing order 4 if 2 isn't
   enough), explicitly set `set_farfield("absorbing")` unless the reference
   design is genuinely a closed box, and check `Z_PROP_BOTTOM`.
4. Run the driven sweep across the same frequency window(s) PalaceForCQED
   used (`run.jl`'s commented-out `driven_lumped(...)` calls at the bottom
   of the file are the likely candidates) and plot `|S11|`/`arg(S11)`
   against `driven_new.svg` directly — visual/numerical dip-depth and
   dip-frequency comparison, not a single scalar like Z0.
5. Document results the same way this session did: what matched, what
   didn't, and for anything that didn't, which single knob (order, boundary,
   geometry dimension, mesh) explained it — verified by an isolated test
   changing only that one thing, not multiple simultaneous changes.
6. Remember the standing rule: **no PALACE/FEM solves locally, ever** — every
   run goes to the CNEA SGE cluster via `qsub` (see project memory
   `local-machine-heavy-fem-jobs-unsafe` and `cluster-frontend-no-direct-work`).
   Check job progress/memory early (`grep "Estimated current\|Estimated peak"`
   on the log) before letting an order-4 run sit unattended — this session
   hit real OOMs from finer-mesh/higher-order combinations more than once.

## Where this session's crosscheck code lives, for reference/reuse

- `impedance_crosscheck/mesh_si_cpw_line.jl`, `gen_mesh.jl`,
  `si_cpw_line.json`, `stage2_analyze.py` — the raw-Gmsh/Julia bare-line
  model and its NRW Z0 extraction.
- `driven_sweep/check_line_impedance.py`, `driven_sweep/submit_line_
  impedance.sh` — the SQDMetal bare-line model, now with `MATCH_JULIA=1`
  mode and env-overridable `SOLVER_ORDER`, `FARFIELD_TYPE`, `AMR_ITERS`/
  `AMR_TOL`, `LEN_LAUNCH_UM`, `GMSH_DIST_SAMPLING`, mesh taper/size knobs —
  useful as a reference for which SQDMetal API calls matter and how they're
  wired, even though the actual task here needs a different (much more
  complex) geometry.
