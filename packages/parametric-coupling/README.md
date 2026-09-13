# parametric-coupling

Hamiltonian-parameter extraction and analytic design helpers for the group's
SQDMetal / Palace transmon-coupler work.

Modules (import as `from parametric_coupling import ...`):

- `DesignHelper` — analytic formulas: transmon EJ/EC/frequency/anharmonicity,
  CPW geometry (impedance, kinetic inductance), resonator Qc, dispersive shifts,
  Purcell rates, and empirical capacitance models fitted to Palace sweeps
  (`circular_island_capacitance`, `interdigital_coupler_capacitance` and their
  inverses).
- `Chip`, `Component` / `Resonator` / `Qubit` / `GroundPlane` / `FloatingSQUID` —
  a small circuit-graph layer that turns capacitance / inductance matrices into
  Hamiltonian parameters.

## Install (editable, into whichever conda env)

```bash
conda run -n <env> pip install -e /home/ruiz/Documents/thesis/packages/parametric-coupling
```

Currently installed editable into `qcg-quantum-design`.
