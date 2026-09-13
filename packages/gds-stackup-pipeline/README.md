# gds-stackup-pipeline

Runs SQDMetal PALACE eigenmode simulations from a raw GDS layout + an XML
material stackup, via SQDMetal's native `gds_design=` input branch. GDS-input
counterpart to the QDesign-based pattern in
`notebooks/design/squid_array_coupling/build_and_eigenmode.py`.

Unlike `packages/parametric-coupling/` and `packages/qcg-quantum-design/`,
this package is **not** simulation-backend-agnostic — it imports SQDMetal
directly, since gluing GDS+stackup input into SQDMetal's API is its whole
job. It only installs into envs that have SQDMetal (py3.11-pinned, e.g.
`qcg-quantum-design`/`qcg-general-design`), not into `quantum`.

See `notebooks/design/gds_stackup_spike/COMPARISON_SQDMETAL_VS_GDS2PALACE.md`
for the source-grounded background (why GDS ingestion needed no new SQDMetal
code, what's actually missing, and the decided port-marking convention:
explicit port geometry passed by the caller, no GDS-layer convention).

## Install

```
conda run -n qcg-general-design pip install -e packages/gds-stackup-pipeline/
```

## Scope

Single-substrate, single-metal-layer stackups only (matches SQDMetal's own
`GeomGDS` ceiling — one chip thickness, one material). Eigenmode simulations
only. Multi-layer stackups and other Palace simulation types (driven-RF,
capacitance) are not supported — extend `stackup.py`/`build.py` if needed.

## Self-check

```
conda run -n qcg-general-design python -m gds_stackup_pipeline.test_stackup
```

## Gotchas (from the first real qpdk cluster run)

Hit while running a genuine qpdk-provided geometry (launcher + coupled
CPW resonator) through this pipeline — none are bugs in this package
itself, but each one silently breaks a run if unaddressed, so they're
recorded here rather than left to be rediscovered.

- **A `gds.write_gds()`'d file from gdsfactory reads as empty.** gdsfactory's
  default `write_gds()` prepends a `$$$CONTEXT_INFO$$$` metadata pseudo-cell
  as `cells[0]`, pushing the real top cell to index 1. `GeomGDS` defaults to
  `gds_cell_index=0` (correct for plain `gdstk`/`gdspy`-written files, which
  don't have this cell) — so it silently reads zero real polygons instead of
  erroring. Fix: `component.write_gds(path, with_metadata=False)`.

- **`create_port_CPW_via_edge_point` fails with "CPW couldn't be formed"**
  (or, worse, silently finds 0 corners with no crash at all until the final
  assert) if the reference point sits *exactly* on the metal boundary rather
  than strictly inside or outside it. The algorithm computes a perpendicular
  direction as `(nearest_boundary_point - reference_point)`, normalized —
  exactly-zero distance makes that a `NaN` vector, and every downstream
  "is this corner collinear" check then silently reads false. Confirmed this
  can even appear to work by accident on one GEOS version (floating-point
  noise avoiding a true zero) and fail on another (GEOS 3.10.6 vs 3.13.1,
  same input, same SQDMetal source) — not a GEOS-version bug to work around,
  the input just wasn't a valid "near the edge, not on it" point to begin
  with. Fix: inset the reference point a few µm into the metal, opposite the
  port's own outward-facing orientation if it's a gdsfactory port.

- **"Could not create line" inside `construct_geometry_in_GMSH`** when the
  GDS's metal bounding box is used as the chip/substrate size verbatim
  (the default when `segment_rectangle` isn't given): every metal edge then
  touches the substrate boundary exactly, which can produce a degenerate
  (near-zero-length) polygon edge in SQDMetal's substrate-minus-metal gap
  computation. A metal-simplification `threshold` large enough to dodge the
  crash is *also* large enough to bridge a real narrow CPW gap shut —
  confirmed directly, it visibly corrupted the meander in the mesh preview.
  Fix: pass `segment_rectangle` a few µm larger than the metal's own bbox on
  every side, so no edge coincides exactly.

- **`fine_mesh_in_rectangle` meshes the *entire* rectangle interior at
  `mesh_min`, not just near the metal.** GMSH's distance field for a filled
  surface region is zero everywhere inside it — taper only happens outside
  the rectangle. Sizing the rectangle to cover most of the chip (e.g. the
  metal's own bbox) produces a uniformly fine mesh with no taper anywhere
  inside the circuit's own footprint, which was slow enough to blow through
  a 1-hour cluster job budget without converging. `fine_mesh_features
  (max_feature_size=...)` is geometry-aware instead (auto-detects the actual
  narrow trace/gap) and produced a ~3x smaller mesh for the same design,
  finishing in ~44 minutes. Prefer it for CPW-style layouts; see
  `run_gds_eigenmode`'s docstring for both options.

- **`pyvista` (imported by SQDMetal's `PVDVTU_Viewer`, itself imported by
  `Eigenmode_Simulation`) resets the matplotlib backend to a Qt one as a
  side effect of initializing VTK**, which crashes on a headless cluster
  node ("requires the 'qt' interactive framework, as 'headless' is
  currently running") — even after an explicit `matplotlib.use("Agg")`
  *before* importing this package, and even with `MPLBACKEND=Agg` set in
  the job's environment. Confirmed by instrumented trace: backend reads
  `"Agg"` right after an early `matplotlib.use()` call, then flips to
  `"QtAgg"` the moment `gds_stackup_pipeline` (→ SQDMetal → `PVDVTU_Viewer`
  → `pyvista`) is imported, with nothing else in between. Fix:
  `matplotlib.use("Agg", force=True)` *after* importing this package, not
  before — see the top of `run_qpdk_launcher_resonator.py` for the working
  order.
