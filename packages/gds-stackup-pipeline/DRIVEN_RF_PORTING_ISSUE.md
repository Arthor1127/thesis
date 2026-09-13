# Driven-RF extension: geometry ingestion works, port creation doesn't fit this topology

## Context

Testing whether this package's `gds_design=` approach extends from eigenmode
(its only current wrapper) to driven-RF, using a genuine reference geometry:
PalaceForCQED's own `2DQv8_eb4_edit.geo` (a real, hardware-validated device,
see `notebooks/design/parametric_coupling_chip/s11_sqdmetal_validation/
HANDOFF.md` for the full cross-check this supports). Repro files are in
`notebooks/design/gds_stackup_spike/driven_hanger_attempt/`:

- `geo_to_gds.py` — parses the `.geo` file's raw `Point(...)` coordinates for
  a subset of its 22 boundary elements (element 0, a bare hanger resonator,
  and elements 4/5, the feedline) and writes them 1:1 to a GDS file
  (`phase2_geo_resonator.gds`) via `gdstk`. No re-derivation of dimensions —
  literal point-for-point conversion, confirmed correct by re-rendering the
  GDS content and comparing to the original `.geo` plot (matches exactly).
- `build_gds_sim.py` — feeds that GDS into `PALACE_Driven_Simulation`
  directly (`gds_design=...`), bypassing this package's `build.py` entirely
  since it only wraps `PALACE_Eigenmode_Simulation` today. Geometry + mesh
  only (`prepare_simulation()`), no `run()`, no cluster job.
- `gds_mesh_overview.png` — the resulting mesh (see below).

## What works

`GeomGDS` + `add_metallic(1)` + `prepare_simulation()` succeed end-to-end on
this real reference geometry. `gds_mesh_overview.png` shows the resonator
meander and feedline (with its taper detail) meshed correctly, matching the
`.geo` file's actual shapes exactly. **GDS ingestion itself is not the
problem** — this confirms the same headline finding
`COMPARISON_SQDMETAL_VS_GDS2PALACE.md` already made for eigenmode also holds
for driven.

## What doesn't work: no GDS-compatible port constructor fits this geometry

`create_port_CPW_via_edge_point((0e-6, 2040e-6), 20e-6, impedance_R=50.0)`
fails:

```
AssertionError: Seems that a CPW couldn't be formed. Check that the
line-segment closest to the point (0.0, 0.00204) is indeed the edge of the
central conductor of the CPW feed.
```

Root cause, not a bug in the constructor: this reference geometry's feedline
is **two disconnected floating conductors** (elements 4 and 5, a combined
20µm-wide gap between them, each 10µm wide) with **no enclosing ground plane
anywhere** — confirmed by accounting for all 22 elements in the full `.geo`
file (resonators, qubit pads, port taper fingers — none is a ground-fill
sheet). `create_port_CPW_via_edge_point` assumes the standard
ground–trace–ground CPW cross-section (finds the nearest edge, then expects a
symmetric ground on the perpendicular far side); this geometry has no
"central conductor with ground on both sides" for it to find, so it fails
deterministically, not from a bad reference point (already tried moving it
around this session — see the isolated-conductor problem, not a placement
problem).

The two GDS-compatible port constructors are exhaustive here
(`COMPARISON_SQDMETAL_VS_GDS2PALACE.md`'s own table — the other 3 need
QDesign pin metadata, `create_port_2_conds` included, so none of those apply
to raw GDS either):

- `create_port_CPW_via_edge_point` — wrong topology assumption, as above.
- `create_waveport_on_boundary` — different port physics entirely (modal,
  not lumped-R), and would need the domain's boundary wall repositioned to
  coincide exactly with the trace's tip for the port to actually couple to
  the two-conductor line instead of exciting open vacuum. Not attempted —
  real design decision (does a wave-port mode-solve on an ungrounded
  two-strip line even converge on something physically meaningful?), not a
  quick fix.

PalaceForCQED's own config (read directly,
`palace_cqed_reference/2DQv9_eb4_driven_lumped.json`) defines its lumped port
as a raw boundary attribute tag (`"Elements":[{"Attributes":[3001],
"Direction":"+Y"}]`) on a specific mesh surface their own `MeshGenerator.jl`
tags directly — a lower-level mechanism than any of SQDMetal's convenience
port finders, which all assume one of a few standard fixed topologies
(CPW-with-ground, QDesign pin pair, or a domain-boundary mode). There's
currently no way in this package (or, as far as this session's reading went,
in `PALACE_Model_Base`'s own port API) to tag an arbitrary custom GDS/edge
region as a lumped port directly, the way raw Palace config or gds2palace's
own GDS-layer-as-port convention (see `COMPARISON_SQDMETAL_VS_GDS2PALACE.md`,
option (a)) does.

## Secondary, smaller issue: fallback mesh is needlessly heavy

Since porting failed, `build_gds_sim.py` fell back to
`fine_mesh_in_rectangle()` over the whole `segment_rectangle` to still test
whether meshing alone succeeds. This reproduces the README's own already-
documented gotcha exactly: it fine-meshed the *entire* rectangle interior
uniformly (268k nodes, 1.6M elements for a fairly small 2-conductor layout),
not just near the metal. Once a port exists, switch to
`fine_mesh_features(max_feature_size=...)` as the README already recommends
for CPW-style layouts — not a new finding, just confirming it applies here
too.

## Possible directions (not decided — picking one is the actual task)

1. Extend `build.py` with a `run_gds_driven` that only supports the
   ground–trace–ground topology `create_port_CPW_via_edge_point` already
   handles, and treat "no ground plane in the GDS" as an unsupported input
   (document the limitation rather than fixing it).
2. Add a raw "tag this GDS region/edge as a lumped port" mechanism —
   equivalent to gds2palace's own layer-convention port marking
   (`COMPARISON_SQDMETAL_VS_GDS2PALACE.md`'s option (a)), which would need
   new code in `PALACE_Model_Base`/`GeomGDS`, not just this wrapper package.
3. Try the wave-port route properly: reposition the relevant domain boundary
   to coincide with the conductors' tip (`set_yBoundary_as_absolute` or
   similar) and confirm PALACE's mode solver actually finds a sane guided
   mode on an ungrounded two-conductor stripline. Untested this session —
   real uncertainty, not just unwritten code.
