# SQDMetal vs gds2palace — source-grounded comparison

Phase 1 of unifying the GDS/KLayout design path with the existing SQDMetal
path. Every claim below comes from reading the installed source, not either
project's docs/README. Locked direction: SQDMetal's Gmsh/Palace backend
stays the single simulation engine; gds2palace is read only for its XML
stackup schema, never wired in as a second engine.

Sources: SQDMetal at `/home/ruiz/SQDmetal/SQDMetal/SQDMetal` (env
`qcg-quantum-design`/`qcg-general-design`, py3.11); gds2palace v0.4.4 at
`/home/ruiz/miniconda3/envs/qcg-gds-design/lib/python3.12/site-packages/gds2palace/`
(env `qcg-gds-design`, py3.12).

| | SQDMetal | gds2palace |
|---|---|---|
| **Geometry ingestion** | `GeomQiskitMetal`/`GeomGDS` both implement `GeomBase` (`Utilities/GeometryProcessors/`). `GeomGDS.process_layers` (`GeomGDS.py:51-101`) reads GDS via `gdstk`. `PALACE_Model_Base` already accepts `gds_design=<path>` as a first-class alternative to `metal_design=<QDesign>` (`Model.py:66-69`, dispatch `Model.py:175-179`) — **GDS ingestion into SQDMetal already works today**, no new adapter code needed for this part. | `util_gds_reader.py:502` `read_gds()`, via `gdspy`, flattens the cell and returns a flat list of layer-tagged polygon objects (`all_polygons_list`). |
| **Stackup / material spec** | Chip thickness: single scalar kwarg `gds_chip_thickness` (`GeomGDS.py:15`, default 500e-6) — no multi-layer/z-position concept for GDS input. Material: a separate string `dielectric_material`, looked up in a hardcoded table (`Utilities/Materials.py:9-60`, e.g. silicon ε=11.45 at `:26-29`); same mechanism for QDesign or GDS input, not read from either. | Full XML stackup schema, parsed with stdlib `ElementTree` (`util_stackup_reader.py`, `read_substrate`/`parse_substrate` at `:1456,1489`): `<Material Name Type Permittivity DielectricLossTangent Conductivity ...>` (`:540-551`), `<Dielectric Name Material Reference Zmin Zmax Thickness>` (`:628-683`), `<Layer Name Layer Type Material Zmin Zmax>` (`:927-974`, Type ∈ CONDUCTOR/VIA/SHEET/DIELECTRIC). Supports multiple z-stacked dielectrics/metals; SQDMetal supports none of that for GDS input. |
| **Port / pin definition** | 5 port constructors on `PALACE_Model_Base_RF`. 4 need QDesign pin metadata and assert `isinstance(self._geom_processor, GeomQiskitMetal)` (`Model.py:1228,1282,1351,1388`) — **don't work with GDS input**. 2 are geometry-only and do work with GDS: `create_waveport_on_boundary` (`Model.py:1479-1518`, plain x/y/z/plane) and `create_port_CPW_via_edge_point` (`Model.py:1400-1443` → `GeomBase.create_CPW_feed_via_point`, `GeomBase.py:39-131`, pure shapely). | A GDS layer is designated a "port layer" via a Python-built `simulation_port` object (`util_simulation_setup.py:33`); polygons on that layer *absent* from the stackup become ports (`util_simulation_setup.py:768-795`). Not parsed from any file — built by the calling script. |
| **Solver types** | Eigenmode (`PALACE_Eigenmode_Simulation`, confirmed by this spike's output `"Problem":{"Type":"Eigenmode"}`), plus Capacitance and Driven-RF simulation classes exist in the same `PALACE/` package. | Palace `Problem.Type` hardcoded `"Driven"` (`util_simulation_setup.py:1884`) + a separate Elmer-based steady-state thermal path (`create_elmer_thermal`, `:1054`). No Eigenmode/Electrostatic/Magnetostatic problem type anywhere in the source. |
| **Qubit-specific postprocessing (EPR/LOM)** | Has EPR-style postprocessing built into the eigenmode workflow (`_EPR_setup` flag in `Model.py:141`, EPR-specific methods elsewhere in `PALACE/`). | Grep across all `.py` files for `EPR|energy_participation|LOM|anharmonic|qubit`: **zero real hits** (only false-positive substring matches like "representation"). No qubit/EPR/LOM postprocessing exists in gds2palace at all — it's exclusively RFIC driven-S-parameter + thermal. |

## Headline finding

The comparison reframes the task: **GDS geometry ingestion into SQDMetal is
not a gap** — `gds_design=<path>` already works, proven end-to-end by the
spike below. The two real gaps are:

1. **Stackup**: `GeomGDS` supports only one substrate material + one uniform
   chip thickness, as plain kwargs — no file-driven, multi-layer stackup
   concept. gds2palace's XML schema is a good source for what a richer
   stackup *would* need to express, but adapting a genuinely multi-layer
   stackup into SQDMetal's single-substrate assumption is unsolved and out
   of scope for this phase.
2. **Ports**: raw GDS has no netlist, so only the 2 geometry-only port
   constructors are usable. See below.

## Port-marking convention — options, not a decision

How should a future GDS-only design mark port locations, given raw GDS has
no pin/netlist metadata? Three options, presented for the user to choose —
not decided here:

- **(a) gds2palace-style layer convention** — dedicate a GDS layer number as
  a "port layer"; polygons on it (with no stackup entry) become ports, same
  convention gds2palace already uses (`util_simulation_setup.py:768-795`).
  *Pro*: matches an established, layer-native convention. *Con*: needs a
  layer-number registry/convention doc; SQDMetal has no built-in reader for
  this, so a translation layer is still needed to turn a marked layer into
  calls to SQDMetal's own port constructors.
- **(b) SQDMetal-native, geometry-only** — no GDS convention at all; ports
  specified as explicit coordinates via `create_waveport_on_boundary` /
  `create_port_CPW_via_edge_point`, passed alongside the GDS file (e.g. in
  the calling script or a sidecar config). *Pro*: zero new code — this is
  what the spike below already uses, and it works today. *Con*: port
  locations live outside the GDS layout, easy to drift out of sync with the
  geometry as a design changes.
- **(c) Hybrid** — a GDS layer as a human-readable marker, resolved against
  a small sidecar JSON/YAML manifest (layer → port type/impedance/direction)
  that a thin adapter reads and turns into calls to (b). *Pro*: keeps ports
  visible in the layout tool while still routing through SQDMetal's existing
  geometry-only constructors. *Con*: another file to keep in sync; more code
  than (a) or (b).

## Spike result

See `test_stackup.xml`, `parse_stackup_to_json.py`, and
`build_sqdmetal_from_stackup.py` in this directory. Result: **PASS**
(`spike_result.json`) — a hand-written XML stackup (1 material, 1
dielectric, 1 metal layer), parsed with gds2palace's own `read_substrate()`
(reused, not reimplemented), resolved into `dielectric_material`/
`gds_chip_thickness` kwargs SQDMetal already accepts, fed into
`PALACE_Eigenmode_Simulation(gds_design=..., ...)` alongside a trivially
generated test GDS, reached `prepare_simulation()` and produced a real
9630-line Gmsh `.msh` mesh plus a valid Palace eigenmode config JSON
(`spike_sim_output/gds_stackup_spike/`).

One real SQDMetal usage requirement surfaced along the way, unrelated to
GDS: `GMSH_Mesh_Builder.build_mesh()` (`GMSH_Mesh_Builder.py:57`) crashes
with an `UnboundLocalError` if zero fine-mesh regions are registered before
`prepare_simulation()`. QDesign-based workflows normally call
`fine_mesh_components()` (QDesign-only) or similar; the GDS-only spike used
the geometry-only `fine_mesh_in_rectangle()` (`Model.py:751-787`) instead —
this applies to any GDS-only SQDMetal usage, not just this spike.

**Two-env constraint**: gds2palace's source uses PEP 701 nested-quote
f-strings (`util_elmer.py:487` in both v0.4.4 and v0.5.0), which require
Python 3.12+. SQDMetal requires `<3.12`. The two packages cannot share one
conda env — confirmed by attempting exactly that (`qcg-general-design`,
cloned from `qcg-quantum-design`, py3.11): `import gds2palace` raises a
`SyntaxError` on that line. The spike therefore runs as two steps in two
envs with a JSON handoff (`resolved_stackup.json`), not a single script.
