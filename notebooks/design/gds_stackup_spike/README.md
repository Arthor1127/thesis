# gds_stackup_spike

Phase 1 spike: prove (or disprove) a GDS + gds2palace-style XML stackup can
reach SQDMetal's existing Gmsh-building code, no changes to SQDMetal
itself. See `COMPARISON_SQDMETAL_VS_GDS2PALACE.md` for the source-grounded
findings and the port-marking convention options.

## Envs

Two conda envs, JSON handoff — gds2palace needs Python 3.12+, SQDMetal
needs <3.12, so they can't share one process (see comparison doc).

- `qcg-gds-design` (py3.12, has gds2palace) — parses the XML stackup.
- `qcg-general-design` (py3.11, clone of `qcg-quantum-design`, has
  SQDMetal) — builds the SQDMetal model and Gmsh mesh. Built once via:
  ```
  conda create -n qcg-general-design --clone qcg-quantum-design
  ```
  See `scripts/claude/sqdmetal_install_guide.md` for how `qcg-quantum-design`
  itself was set up.

## Run

```
conda run -n qcg-gds-design python parse_stackup_to_json.py
conda run -n qcg-general-design python build_sqdmetal_from_stackup.py
```

Check `spike_result.json` for `"status": "pass"` and
`spike_sim_output/gds_stackup_spike/` for a real `.msh` mesh + Palace
config JSON.
