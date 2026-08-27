# Plan: Reimplement QDesignOptimizer on SQDMetal/PALACE (new package `QCGOptimizer`)

**Naming**: "QCGOptimizer" stands for **Quantum Circuits Group Optimizer**. Repo/directory: `/home/ruiz/QCGOptimizer`. Python import package: `qcgoptimizer` (lowercase, per PEP 8). Internal class names keep the `Palace` suffix (`DesignAnalysisPalace`, `MiniStudyPalace`, `CapacitanceMatrixStudyPalace`, ...) to document that these are the PALACE-backed implementations — only the project/package name changed, not the backend-labeling convention inside the code.

## Context

`qdesignoptimizer` (installed in conda env `qdesignoptimizer`, source at github.com/202Q-lab/QDesignOptimizer, pinned commit `d4f6ada5ada59b786df1006d53f8f148b364364e`) automates quantum-chip design optimization on top of qiskit-metal, but its simulation layer is hardcoded to Ansys: `design_analysis.py` builds `EPRanalysis(design, "hfss")`/`LOManalysis(design, "q3d")`, drives the Ansys COM session directly via a custom `_ComGeometryModeler`, and calls `pyEPR` for nonlinear Hamiltonian (EPR/Kerr/chi) extraction. This is not behind a swappable interface — it's Ansys API calls woven through the orchestration logic.

The user has SQDMetal (env `qcg-quantum-design`) available, which drives the open-source PALACE FEM solver instead of Ansys, and — critically — already depends on the exact same forked qiskit-metal package (`quantum-metal`) that `qdesignoptimizer` itself requires (`Requires-Dist: quantum-metal>=0.7.4`). SQDMetal is not a qiskit-metal renderer plugin like the Ansys HFSS/Q3D renderers; it's a standalone class (`PALACE_Eigenmode_Simulation`, `PALACE_Capacitance_Simulation`, etc.) that takes a live `QDesign` object directly and builds its own Shapely→Gmsh geometry/mesh pipeline, independent of qiskit-metal's renderer system. This makes it a clean, if non-trivial, integration point: qdesignoptimizer's qiskit-metal geometry-building code (chip creation, component placement, `design.variables`, the `render_qiskit_metal(design, **kwargs)` callback) is reusable unchanged; only the simulation-orchestration layer needs replacing.

The goal is a new package, **`qcgoptimizer`**, that reimplements all of qdesignoptimizer's user-facing functionality (target/optimization DSL, iterative optimization loop, eigenmode+EPR extraction, capacitance-matrix decay studies, surface-participation/TLS-loss analysis, plotting) against SQDMetal/PALACE instead of Ansys — without modifying either existing conda environment (`qdesignoptimizer`, `qcg-quantum-design`).

This plan targets all 4 phases (full feature parity, not just the eigenmode core), and builds the new environment as a clone of `qcg-quantum-design` (non-destructive) rather than from scratch.

## Architecture decision: vendor + new orchestrator (not a fork, not a full-package dependency)

Three options were considered:

- **Fork qdesignoptimizer and subclass/monkeypatch** — rejected: `DesignAnalysis.__init__` does all Ansys setup directly in the constructor body with no override seam; a subclass would have to skip `super().__init__()` entirely, buying nothing over a fresh class while adding brittle coupling to upstream internals. It would also leave ~40% of `design_analysis.py` and all of `sim_capacitance_matrix.py` as dead Ansys code in the fork, and collide on the `qdesignoptimizer` import name if both backends are ever needed side-by-side to cross-validate.
- **Depend on `qdesignoptimizer` as a normal pip/git dependency** — confirmed technically safe (its modules import cleanly with no Ansys installed — only `pyaedt`/`pyEPR` touch COM at *call* time, not import time) but rejected as primary: pip would still install the full heavy dependency closet (`pyaedt`, `pyEPR-quantum`, `pyside6`, `qutip`, `scqubits`, `qdarkstyle`) regardless of which submodule is imported, risking version conflicts with SQDMetal's own pins.
- **Vendor the backend-agnostic files verbatim, write a new orchestrator — RECOMMENDED.** Lean, SQDMetal-scoped dependencies; an explicit, auditable reuse boundary; no namespace collision (both packages installable side by side to compare backends); clean room for new orchestration logic with no inherited Ansys-era constructor behavior.

Implementation will need `git clone https://github.com/202Q-lab/QDesignOptimizer && git checkout d4f6ada5ada59b786df1006d53f8f148b364364e` to get raw source files to vendor (the installed package alone doesn't include example/test code that may clarify `jj_setup` shapes — see Risk 3 below).

## Package structure

```
qcgoptimizer/
├── pyproject.toml
├── VENDORED_FROM.md                     # upstream repo + commit + file list + date
├── src/qcgoptimizer/
│   ├── __init__.py
│   ├── _vendor/                         # verbatim copies from qdesignoptimizer@d4f6ada
│   │   ├── anmod_optimizer.py           #   ANModOptimizer — unchanged
│   │   ├── logger.py                    #   unchanged
│   │   ├── sim_plot_progress.py         #   unchanged
│   │   ├── estimation/purcell_limit.py                          # unchanged
│   │   ├── estimation/classical_model_decay_into_charge_line.py # unchanged
│   │   └── utils/{names_parameters,names_design_variables,names_qiskit_components,
│   │              optimization_targets,plotting,chip_generation,utils}.py
│   │              # unchanged, except utils.py loses close_ansys()
│   ├── design_analysis_types.py         # OptTarget/DesignAnalysisState/SimulationResults/
│   │                                     # SurfaceProperties/Interfaces/InterfaceProperties,
│   │                                     # re-exported near-verbatim (MiniStudy/MeshingMap NOT copied)
│   ├── mini_study.py                    # NEW: MiniStudyPalace, FineMeshSpec
│   ├── materials.py                     # NEW: SurfaceProperties -> SQDMetal MaterialInterface mapping
│   ├── sim_capacitance_matrix.py        # NEW: CapacitanceMatrixStudyPalace, ModeDecayStudyPalace (ABC),
│   │                                     #      ModeDecayIntoChargeLineStudyPalace, ResonatorDecayIntoWaveguideStudyPalace
│   ├── design_analysis.py               # NEW: DesignAnalysisPalace orchestrator;
│   │                                     #      save_optimization_results/merge_partitioned_simulation vendored verbatim
│   ├── epr_mapping.py                   # NEW: pure functions, SQDMetal EPR/eig output -> system_optimized_params
│   └── port_mapping.py                  # NEW: pure functions, jj_setup/port_list -> SQDMetal create_port_* calls
├── tests/
│   ├── test_epr_mapping.py              # synthetic DataFrames, no Palace binary needed
│   ├── test_port_mapping.py
│   ├── test_capacitance_indexing.py
│   └── integration/test_transmon_res_eigenmode.py   # needs real Palace binary; skip-marked
├── examples/transmon_resonator/         # mirrors SQDMetal's Palace_Eigenmode_Transmon_Res.ipynb
└── docs/
```

Each `_vendor/` file carries a header noting origin commit + license (verify Apache-2.0 at clone time), so future re-vendoring against a newer upstream commit is a simple diff.

## MiniStudyPalace — replacing MiniStudy field-by-field

New dataclass in `mini_study.py`. Kept identical (name + semantics) to upstream `MiniStudy` so `OptTarget`/`system_target_params` user code needs no changes: `qiskit_component_names`, `port_list`, `open_pins`, `modes`, `jj_setup`, `design_name`, `adjustment_rate`, `render_qiskit_metal_eigenmode_kw_args`, `run_capacitance_studies_only`, `capacitance_matrix_studies`, `surface_properties`.

Replaced/new fields:

| Old HFSS field | New PALACE field | Maps to |
|---|---|---|
| `project_name` | dropped | n/a |
| `x_buffer_width_mm`, `y_buffer_width_mm` | `boundary_proportion_xy`, `boundary_proportion_z` | `PALACE_Model_Base.set_xBoundary_as_proportion` / `set_yBoundary_as_proportion` / `set_zBoundary_as_proportion` |
| `max_mesh_length_port`, `max_mesh_length_lines_to_ports`, `hfss_wire_bond_*`, `build_fine_mesh` | `fine_mesh_specs: List[FineMeshSpec]` | `fine_mesh_components()` / `fine_mesh_along_path()` / `fine_mesh_in_rectangle()` |
| `nbr_passes` | dropped (warn if set) | no adaptive-pass concept in a direct FEM solve |
| `delta_f` | `solver_tol` | `palace_user_options["solver_tol"]` |
| `cos_trunc`, `fock_trunc` | dropped | SQDMetal's EPR path is a perturbative closed-form (Chi from EPR+Ej), no Fock-truncation knob |
| — | `palace_user_options: dict` | passed straight through to `PALACE_*_Simulation(user_options=...)`: `starting_freq`, `number_of_freqs`, `solns_to_save`, `solver_order`, `solver_maxits`, `dielectric_material`, `mesh_min/max`, `taper_dist_min/max`, `fillet_resolution`, `palace_mode`, `palace_dir`, `num_cpus`, `num_threads`, etc. |
| — | `mode: "PC"\|"simPC"\|"HPC"` | `PALACE_*_Simulation(mode=...)` |
| — | `meshing: "GMSH"\|"COMSOL"` | `PALACE_*_Simulation(meshing=...)` — GMSH default, no license needed |
| — | `sim_parent_directory: str` | where Palace input/output land |
| — | `epr_interface_thicknesses: dict` | passed to `setup_EPR_interfaces()` |

`FineMeshSpec` (`kind: "components"|"path"|"rectangle"` + `kwargs` bag) replaces `MeshingMap` entirely — SQDMetal's fine-mesh API is region/component-driven, with no per-component-class name-generator equivalent.

## Eigenmode + EPR: `run_eigenmodes()` / `run_epr()`

Mirrors SQDMetal's `Palace_Eigenmode_Transmon_Res.ipynb` call order, in `DesignAnalysisPalace`:

1. `self.render_qiskit_metal(self.design, **mini_study.render_qiskit_metal_eigenmode_kw_args)` — unchanged, no `_ComGeometryModeler` needed at all.
2. `eigen_sim = PALACE_Eigenmode_Simulation(name=..., sim_parent_directory=..., mode=..., meshing=..., user_options=palace_user_options, metal_design=self.design)`.
3. `eigen_sim.add_metallic(1); eigen_sim.add_ground_plane()`.
4. `port_mapping.add_jj_ports(eigen_sim, mini_study.jj_setup)` → `eigen_sim.create_port_JosephsonJunction(qObjName, L_J=..., C_J=...)`. **Real translation, not pass-through**: upstream's `jj_setup` names raw HFSS rect/line objects; `create_port_JosephsonJunction` instead reads `design.qgeometry.tables['junction']` by component id + `junction_index`. Resolve exact upstream `jj_setup` shapes from the cloned repo's examples/tests in Phase 1 (see Risk 3).
5. `port_mapping.add_ports_from_port_list(eigen_sim, mini_study.port_list, design)` — dispatches per qiskit-metal component type to `create_port_CPW_on_Launcher`/`create_port_CPW_on_Route` with the given impedance; this is what makes kappa extractable from the port-EPR output.
6. If `mini_study.surface_properties`: `eigen_sim.setup_EPR_interfaces(...)` (see below).
7. Apply `mini_study.fine_mesh_specs` via the matching `eigen_sim.fine_mesh_*` calls.
8. `eigen_sim.prepare_simulation(); eigen_sim.run()` (blocking local subprocess).
9. Retrieve: `eigen_sim.retrieve_data()` (freq/Q), `eigen_sim.retrieve_mode_port_EPR()` (`mat_mode_port`, `eigenfrequencies`, `loaded_Q`, `kappa`), and if junctions exist, `eigen_sim.calculate_hamiltonian_parameters_EPR(print_output=False)` (`f_modes_GHz`, `f_norms_GHz`, `EPR`, `Chi`, `Lamb`, `Detuning`).

`run_eigenmodes()` and `run_epr()` stay as two public methods for API-shape parity, but internally `run_epr()` reads back data already produced by the single Palace call `run_eigenmodes()` triggered (cached as `self._eigen_sim`) — merging the two Ansys-era simulation calls into one Palace run.

**Mapping to `system_optimized_params`** (pure functions in `epr_mapping.py`, unit-tested against synthetic DataFrames):
- `freq`: `f_norms_GHz*1e9` (dressed, Lamb-shift-corrected) when EPR ran; else raw `eigenfrequencies` — mirrors upstream's use of pyEPR's dressed frequencies.
- `kappa`: from `retrieve_mode_port_EPR()['kappa']` (native `port-Q.csv` byproduct of the same eigenmode run), **not** a separate driven S-parameter sweep — keeps the "one simulation per outer iteration" cost model intact and is the direct physical analogue of upstream's `Freq/Q` kappa. A driven-sweep fallback (`PALACE_Driven_Simulation`) is a documented future option, not the default.
- `nonlinearity`: directly from `hamiltonian['Chi']` (diagonal = anharmonicity, off-diagonal = cross-Kerr, matching upstream's convention) — convert MHz → Hz (`*1e6`).
- `charge_line_limited_t1`: unaffected, comes entirely from the capacitance-matrix path (below).
- Mode↔index alignment: reuse upstream's pattern — sort `mini_study.modes` by target frequency ascending, assume Palace returns modes in ascending-frequency order (true by construction). **Assert this monotonicity explicitly in code** — this is Risk 1 below.
- Junction participation (`param_participation_ratio`): directly from `mat_mode_port` — a structural improvement over upstream (no custom geometry rendering needed).

## Capacitance-matrix studies

New `sim_capacitance_matrix.py`, mirroring upstream's class hierarchy (`CapacitanceMatrixStudyPalace`, `ModeDecayStudyPalace` ABC, `ModeDecayIntoChargeLineStudyPalace`, `ResonatorDecayIntoWaveguideStudyPalace`) so the vendored decay-rate math (`estimation/classical_model_decay_into_charge_line.py`, unchanged) needs zero edits. Only `simulate_capacitance_matrix()` differs: build `PALACE_Capacitance_Simulation(metal_design=design, ...)`, `add_metallic`/`add_ground_plane`, apply fine-mesh specs, `prepare_simulation()`, `run()` → raw NumPy matrix (Farads).

**Row/col indexing risk (flagged, real work item)**: Palace's `terminal-C.csv` indexes conductors by integer position (`Terminal.Index`), with no semantic name — unlike upstream's Q3D/LOManalysis matrix, which is already name-indexed. Handling:
- Always call `cap_sim.display_conductor_indices(save=True)` and persist the figure as run provenance (not optional debugging).
- Build a deterministic name→index adapter (`_infer_conductor_names(cap_sim, design)`, genuinely new code) by matching each contiguous-metal-mapping Shapely polygon against qiskit-metal component/pin footprints from `design.qgeometry.tables`.
- Wrap the raw matrix in a `pd.DataFrame` indexed by inferred names so vendored `.loc[...]`-based downstream code (both `_update_optimized_params_capacitance_simulation` and the decay-study formulas) works unmodified.
- For canonical floating-transmon(+resonator/+feedline) layouts, prefer SQDMetal's own built-in `calc_params_floating_Transmon(_from_files)`/`calc_params_2_floating_transmons(_from_files)` extractors (fixed positional `conductor_indices`) instead of the general name-inference adapter — lower risk for the common case.
- Preserve upstream's defensive `np.abs(...)` around every matrix access (sidesteps most sign-convention risk between Q3D's and Palace's capacitance conventions).

## Surface participation / TLS loss

**Reimplement using SQDMetal's native `setup_EPR_interfaces()`/`retrieve_interface_EPR_data()`**, not a port of upstream's raw Ansys-COM `_surface_rendering_for_surface_participation_ratios()` — Palace computes this as a built-in postprocessing feature (`Boundaries.Postprocessing.Dielectric`), no custom geometry code needed.

- `SurfaceProperties`/`Interfaces`/`InterfaceProperties` dataclasses reused unchanged. Map `interfaces.substrate_air`→`setup_EPR_interfaces(substrate_air=...)`, `interfaces.metal_substrate`→`substrate_metal=...` (argument-order rename only), `interfaces.metal_air`→`metal_air=...`. `interfaces.underside_air` has **no SQDMetal equivalent** — drop/warn, document as a known gap.
- New `materials.py::to_material_interface(props: InterfaceProperties) -> MaterialInterface` converts `eps_r`/`tan_delta_surf`/`th` (mm→m).
- `SurfaceProperties.sheet_material` (Ansys-material-library comment dropped) is inert in the PALACE path — Palace treats metal as ideal PEC by default; finite sheet impedance is a separate opt-in (`add_kinetic_inductance()`, unrelated to named materials). Keep field for forward-compat, document as currently unused.
- `DesignAnalysisPalace.get_surface_p_ratio()` wraps `retrieve_interface_EPR_data()`'s per-mode `{'SA':{p,Q},'MS':{...},'MA':{...}}` into the same `{interface: {mode_idx: value}}` shape upstream returns, so vendored `get_surface_p_ratio_df()` keeps working. Junction participation comes from `mat_mode_port` (already covered above). Upstream's bulk-dielectric participation (`_get_dielectric_p_ratio`) has no direct SQDMetal equivalent (folded into material `LossTan` instead) — resolve as an open design question in Phase 3 (approximate via `1 - Σ(surface participations)`, or drop and document).
- Upstream's `hfss_wire_bonds must be False for surface-p-ratio` assertion has no analogue — wire bonds/air bridges are out of scope entirely (below).

## Explicitly dropped / deferred

- `_ComGeometryModeler` — no replacement needed; SQDMetal builds geometry itself from QGeometry tables.
- `close_ansys()` — dropped from vendored `utils.py`; not needed for a blocking subprocess call.
- HFSS adaptive-pass knobs (`nbr_passes`, `delta_f`, `update_nbr_passes*`, `update_delta_f`) — dropped; nearest analogue (`solver_tol`) folded into `palace_user_options`. Palace's `Model.Refinement.UniformLevels` (uniform, not adaptive) could be exposed later but isn't a semantic equivalent.
- Air bridges / wire bonds (`get_air_bridge_coordinates`, `hfss_wire_bond_*`, COM `create_bondwire`) — dropped, no SQDMetal equivalent found; document as a known limitation if ever needed.
- `cos_trunc`/`fock_trunc` and strongly-nonlinear (e.g. fluxonium) junction treatment — explicitly out of scope, matching SQDMetal's own documented limitation (perturbative EPR only; full numerical diagonalization, e.g. via `scqubits`, would be a separate manual step, not part of this package).
- `project_name` — dropped, no replacement.
- HPC job submission/monitoring — `mode="HPC"` on SQDMetal classes only generates a `.sbatch` file; this plan targets `mode="PC"` (local blocking run) for all phases. Cluster orchestration is a distinct, separately-scoped follow-up if ever needed.
- COMSOL meshing path — deferred; `meshing="GMSH"` only throughout this plan.

## Environment plan

Create a new conda env `qcgoptimizer` as a **clone** of `qcg-quantum-design` (non-destructive — leaves both existing envs untouched):

```
conda create --name qcgoptimizer --clone qcg-quantum-design
conda activate qcgoptimizer
pip install -e /home/ruiz/SQDmetal/SQDMetal        # re-affirm editable SQDMetal install in the new env
pip install -e <path-to-qcgoptimizer-checkout>        # once qcgoptimizer's pyproject.toml exists
```

This inherits SQDMetal/quantum-metal/gmsh/pandas/shapely already known to work together in `qcg-quantum-design`, avoiding a from-scratch dependency resolution. `qcgoptimizer`'s own added dependencies (`pandas`, `numpy`, `matplotlib`, `scipy` for `ANModOptimizer`'s `scipy.optimize.minimize`) are all already satisfiable from the clone.

**Update (done)**: the `qcgoptimizer` env was created via `conda create --name qcgoptimizer --clone qcg-quantum-design`, then `pip install -e /home/ruiz/SQDmetal/SQDMetal --no-deps` and `pip install -e /home/ruiz/QCGOptimizer --no-deps`. Both drift items flagged below turned out to be non-issues in practice:
- `quantum-metal` is `0.7.4` in both `qcg-quantum-design` and the new `qcgoptimizer` clone (the earlier `0.8.1` reading was a research-agent error, not real drift) — matches what `qdesignoptimizer` and SQDMetal both target, no version pin needed.
- `numpy 2.4.6` is present in both envs despite SQDMetal's `constraints.txt` nominally pinning `numpy<2`, but the user's own existing SQDMetal-based scripts (`scripts/python/` in the thesis repo) already run successfully against this exact environment, so this is not a real blocker.
- Verified: `qcgoptimizer`, `SQDMetal.PALACE.Eigenmode_Simulation`, and `SQDMetal.PALACE.Capacitance_Simulation` all import together cleanly in the new env; the original `qcg-quantum-design` env was re-checked afterward and is unaffected (still `quantum-metal==0.7.4`, `sqdmetal==0.0.1.dev1` installed from the same shared editable source).

**Blocker to flag explicitly**: SQDMetal invokes Palace as an external MPI binary via subprocess (`palace_dir` + `palace_mode='local'`); no `palace` binary is currently found on this machine (`which palace` empty). Building/installing it (typically via Spack, per SQDMetal's own HPC docs) is a prerequisite outside this plan's scope but is a hard blocker for any real `run()` call. Phase 1's geometry/config-generation steps (through `prepare_simulation()`) can be validated without it; only the actual solve needs the binary.

## Phased build order

**Phase 1 — Vendor + single-transmon+resonator eigenmode-only loop** (no capacitance, no surface EPR).
Clone upstream at the pinned commit, populate `_vendor/`, write `VENDORED_FROM.md`. Implement `design_analysis_types.py`, `mini_study.py`, `port_mapping.py`, `epr_mapping.py`, and a minimal `DesignAnalysisPalace` (eigenmode+EPR only, `capacitance_matrix_studies=[]`). Build a design mirroring `Palace_Eigenmode_Transmon_Res.ipynb` (transmon + readout resonator, one JJ port, one 50Ω CPW-route port) with targets: qubit freq (via Lj), qubit anharmonicity, resonator freq (via length), resonator kappa (via coupling geometry). Goes first because it exercises the highest-uncertainty new code (EPR mapping, port creation, mode-index alignment) with no compounding risk from capacitance indexing or EPR-interface setup.

**Phase 2 — Capacitance-matrix studies + charge-line/waveguide decay.**
Implement `sim_capacitance_matrix.py` and the conductor-name-inference adapter. Extend the Phase-1 design with a charge line; cross-check against SQDMetal's own `calc_params_floating_Transmon` where the design fits a canonical 3–5-conductor shape.

**Phase 3 — Surface-participation / TLS-loss EPR.**
Implement `materials.py`, wire `setup_EPR_interfaces()` into `run_eigenmodes()`, implement `get_surface_p_ratio()`/`get_surface_p_ratio_df()`. Resolve the bulk-dielectric-participation open question.

**Phase 4 — Polish.**
Verify `save_optimization_results`/`merge_partitioned_simulation` wiring (should be pure copy-paste), verify `plot_progress()`'s expected input shape still matches, docs/examples, partitioned-optimization example if time permits.

## Highest risk / most uncertain steps

1. **Mode-index alignment** between `mini_study.modes` (sorted by target freq) and Palace's eigenvalue-solver mode ordering — silent-failure risk (wrong mode gets the wrong freq/chi). Mitigation: assert monotonicity of returned `eigenfrequencies`, cross-check nearest-simulated-frequency assignment against the naive sorted-index assignment and warn on disagreement.
2. **Capacitance-matrix conductor-index → semantic-name inference** — genuinely new code, no upstream pattern to reuse. Budget the most implementation time here; lean on SQDMetal's fixed-index extractors for canonical designs as a lower-risk fallback.
3. **`jj_setup` dict-shape translation** to `create_port_JosephsonJunction`'s `qObjName`+`junction_index` model — upstream's format is keyed around raw HFSS rect/line names with no meaning in the qiskit-metal junction table. Must validate against the cloned repo's actual example/test designs during Phase 1, since the installed package alone doesn't include them.
   **Update (resolved)**: cloned upstream at the pinned commit and read `src/qdesignoptimizer/utils/names_design_variables.py::junction_setup()`, the function every tutorial (`tutorials/examples_coupled_transmon_chip/mini_studies.py` etc.) uses to populate `jj_setup`. Real shape:
   ```python
   jj_name = f"jj_{name_mode(mode)}"
   setup = {jj_name: {
       "rect": f"JJ_rect_Lj_{name_mode(mode)}_rect_jj",   # Ansys-only, not needed by SQDMetal
       "line": f"JJ_Lj_{name_mode(mode)}_rect_jj_",        # Ansys-only, not needed by SQDMetal
       "Lj_variable": design_var_lj(mode),                 # design.variables key holding L_J
       "Cj_variable": design_var_cj(mode),                 # design.variables key holding C_J
   }}
   ```
   So `port_mapping.add_jj_ports()` doesn't need to parse the `rect`/`line` strings at all — it only needs `Lj_variable`/`Cj_variable` (look up their current values in `design.variables`) plus the qiskit-metal **component name that owns the junction**, which is `name_mode(mode)` itself (the same string `jj_setup`'s key is derived from) — i.e. the mode name already *is* the qiskit-metal component instance name by this codebase's own naming convention, so `create_port_JosephsonJunction(qObjName=name_mode(mode), junction_index=0, L_J=..., C_J=...)` is a direct, low-risk translation. This downgrades what was flagged as the plan's highest-uncertainty item to a straightforward adapter.

## Verification strategy

**Phase 1**: Unit-test `epr_mapping.py` against synthetic `retrieve_mode_port_EPR`/`calculate_hamiltonian_parameters_EPR` return shapes (from documented docstrings) — assert correct keys/units (Hz, not GHz/MHz). Geometry-only smoke test: render → ports → `prepare_simulation()` → confirm a valid mesh + Palace JSON config is produced (`json.load()` succeeds) — runs today, no Palace binary needed. End-to-end (once the binary exists): run several `ANModOptimizer` outer iterations and confirm qubit frequency converges toward target within bounds; cross-check absolute numbers against SQDMetal's own tutorial order-of-magnitude values and closed-form transmon formulas (`f_q ≈ √(8·E_J·E_C)/h − E_C/h`), since no Ansys ground truth exists on this machine.

**Phase 2**: Regression-test the name-inference adapter against a fixture design with manually-verified conductor mapping. Cross-check `ModeDecayIntoChargeLineStudyPalace`/`ResonatorDecayIntoWaveguideStudyPalace` outputs against SQDMetal's own `calc_params_floating_Transmon` fields for the same design (order-of-magnitude/scaling agreement expected, not exact — different derivations).

**Phase 3**: Sanity-check `get_surface_p_ratio_df()` relative magnitudes (junction participation should dominate for a typical transmon) against qualitative expectations from TLS-loss/EPR-method literature — plausibility check, not numeric match.

**Phase 4**: Regenerate `plot_progress()` for a full run and visually confirm convergence. Run `merge_partitioned_simulation` end-to-end and confirm merged output structure matches a single-partition run's (structural regression test).

## Critical files referenced

- `/home/ruiz/miniconda3/envs/qdesignoptimizer/lib/python3.11/site-packages/qdesignoptimizer/anmod_optimizer.py`
- `/home/ruiz/miniconda3/envs/qdesignoptimizer/lib/python3.11/site-packages/qdesignoptimizer/design_analysis.py`
- `/home/ruiz/miniconda3/envs/qdesignoptimizer/lib/python3.11/site-packages/qdesignoptimizer/design_analysis_types.py`
- `/home/ruiz/miniconda3/envs/qdesignoptimizer/lib/python3.11/site-packages/qdesignoptimizer/sim_capacitance_matrix.py`
- `/home/ruiz/SQDmetal/SQDMetal/SQDMetal/PALACE/Eigenmode_Simulation.py`
- `/home/ruiz/SQDmetal/SQDMetal/SQDMetal/PALACE/Capacitance_Simulation.py`
- `/home/ruiz/SQDmetal/SQDMetal/SQDMetal/PALACE/Model.py`
