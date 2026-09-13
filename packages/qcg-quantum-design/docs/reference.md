# `qcg-quantum-design` — reference

Second-quantized model of an **SSH chain of transmons** (capacitive +
cavity-mediated inductive coupling) sharing a **microwave cavity**, plus the
Bogoliubov / normal-ordering engine that diagonalizes it and the tooling to read
physical quantities off the result.

Importable modules (flat, editable-installed — see
[`README.md`](../README.md)):

| module | contents |
|---|---|
| `ssh_chain_sq` | everything below |
| `sq_engine` | shim re-exporting the engine half (`TensorTerm`, `Superposition`, `QuantumOperator`, `SecondQuantizationEngine`, `CompositeEngine`) |
| `thesis_cmaps` | registers the `thesis_cmap` (sequential) and `thesis_cmap_diverging` (blue–white–red) matplotlib colormaps, each with a `_r` reverse; import for side effect, or use `thesis_cmap`, `thesis_cmap_diverging`, `THESIS_PALETTE`. `qcg` / `qcg_r` remain as legacy aliases |
| `thesis_style` | categorical colour cycle (`CYCLE`, `COLORS`, `thesis_cat`), standalone rcParams style `use(typography=…, bw=…)`, and helpers `figsize`, `finalize`, `panel_labels`, `highlight`, `save` |
| `qcg_cmaps` | legacy shim re-exporting `thesis_cmaps` |

```python
from ssh_chain_sq import PHOTONICChain, SecondQuantizationEngine
```

---

## 1. Physical model

`n = 2N` transmon sites, each a bosonic mode `b_j` (or a spinless fermion — see
`statistics`), plus one cavity mode `a`. All energies are in units of the bare
site (qubit) frequency, which is therefore `1`.

**Capacitive (SSH) sector** — number–*non*-conserving quadratic form coming from
the `½CV²` charging energy in the phase/charge basis:

```
H_C = Σ_ij  A_ij b†_i b_j  +  Σ_ij B_ij b†_i b†_j  + h.c.
A = ½(I + T),   B = -½ T
T = alternating nearest-neighbour hopping,  gc1 / gc2 = ¼ Ω_C (1 ± χ_C)
```

`Ω_C` sets the capacitive inter-site coupling scale, `χ_C` the SSH dimerization
(bond staggering). `PBC=True` adds the wrap-around `gc2` bond.

**Josephson / inductive sector** — the cavity-mediated inductive coupling,
`gj1 / gj2 = ± Ω_J (1 ± χ_J)`, assembled into Nambu blocks `AJ, BJ`. In the full
Hamiltonian this sector multiplies the cavity operator `sin(γ(a+a†))`.

**On-site Kerr** — transmon anharmonicity `−kerr Σ_i b†_i b†_i b_i b_i`
(bosons only).

**Cavity** — mode `a` of frequency `omega_r`, truncated to `fock_photon` Fock
levels, coupled through `γ` (the phase zero-point fluctuation, i.e. the argument
of the `sin`).

Full Hamiltonian actually assembled by the class:

```
H = omega_r a†a ⊗ I_matter
  + I_photon ⊗ H_C_fock
  − sin_op(a) ⊗ H_J_fock
```

restricted to a fixed matter-excitation sector, with
`sin_op = γ(a+a†)` (`linear_coupling=True`) or the exact
`sin(γ(a+a†))` (`linear_coupling=False`).

---

## 2. Algebraic data structures

### `TensorTerm(tensor, daggers)`
A single operator string:
`Σ_{i₁…i_k} tensor[i₁,…,i_k] · O_{i₁}^{d₁} … O_{i_k}^{d_k}`, where each
`daggers[m] ∈ {0,1}` marks index `m` as an annihilation (`0`) or creation (`1`)
operator. `tensor` is a rank-`k` `np.ndarray`; `daggers` a length-`k` tuple.
E.g. a Kerr term is `TensorTerm(V, (1,1,0,0))` ↔ `Σ V_{ijkl} b†_i b†_j b_k b_l`.

### `QuantumOperator(terms=None)`
A list of `TensorTerm`s, i.e. a sum of operator strings of possibly different
rank. Operators:

| op | meaning |
|---|---|
| `A + B` | concatenate term lists |
| `scalar * A` | scale every term's tensor |
| `A * B` | operator product via tensor outer product (`rank_A + rank_B`) |
| `A.simplify(tol)` | merge terms with identical `daggers` (add tensors), drop `< tol` |

Algebra is *formal* — no commutators are applied until `normal_order`.

### `Superposition(state_dict=None)`
A quantum state as `{fock_occupation_tuple: complex_amplitude}`. Supports `+`,
left scalar `*`, `.normalize()`, and `.trim(threshold)` (keep the largest-weight
Fock states until cumulative probability ≥ `threshold`, then renormalize).

---

## 3. `SecondQuantizationEngine`

```python
SecondQuantizationEngine(n_modes, statistics='boson',
                         coeff_tol=1e-9, zero_tol=1e-12, use_gpu=False)
```
Fixed number of modes, bosonic or fermionic (Jordan–Wigner phases handled
internally). `use_gpu=True` uses CuPy for the squeezed-vacuum expansion when
available. **Scope:** basis-independent Fock-space algebra + Bogoliubov theory;
it knows nothing about the SSH chain.

### Operator ↔ state

| method | what it does / physical scope |
|---|---|
| `apply_operator(op, superposition)` | act a `QuantumOperator` on a `Superposition`, returns the new `Superposition`. |
| `normal_order(op)` | move all creation operators left using (anti)commutators (Wick), folding rank-`k` → rank-`k−2` → … Returns a normal-ordered `QuantumOperator`. This is where operator identities actually get used. |
| `change_basis(op, U, V)` | substitute `b_i = Σ_j U_ij b_j + V*_ij b†_j` (a Bogoliubov transform) into every term and re-expand. Rank is preserved, term count grows by `2^k`. |
| `exp_operator(op, order)` | Taylor-expand `exp(op)` to `order`, normal-ordered. |

### Bases and matrices

| method | what it does |
|---|---|
| `generate_basis(n_excitations)` | all Fock states with exactly `n_excitations` quanta (stars-and-bars for bosons, `C(n_modes, k)` for fermions). |
| `build_operator_matrix(op, basis)` | dense matrix of `op` in `basis` (rows = bra). |
| `build_sparse_operator_matrix(op, basis)` | same, CSR; bosons use a vectorized batch-apply fast path. |
| `build_quadratic_sparse_matrix(op, basis)` | specialized fast builder for purely quadratic `op`. |
| `build_rectangular_operator_matrix(op, ket_basis, bra_basis)` | matrix of an operator that maps between two *different* excitation sectors (e.g. a pair-creation term). |
| `build_fast_sparse_matrix(op, basis)` | generic COO builder (simplifies first). |
| `build_basis_change_matrix(source_states, expand_fn)` | assemble the (generally rectangular) matrix whose columns are `expand_fn(state)` for each source state; row basis discovered as the union of all expansions. |

### Quadratic diagonalization

| method | what it does / scope |
|---|---|
| `diagonalize_quadratic(op)` | Bogoliubov-diagonalize a quadratic `op` (must contain only `(1,0)` and `(1,1)` terms). Bosons: Van Hemmen / Cholesky + symplectic Gram–Schmidt, positive-frequency branch. Fermions: Hermitian Nambu matrix, particle-block selection. Returns `(frequencies, U, V)` with `U,V` shape `(n_modes, n_modes)`. |
| `build_squeezed_vacuum(U, V, truncation, threshold=0.99)` | the analytic vacuum annihilated by the Bogoliubov modes, `exp(−½ Σ M_ij b†_i b†_j)|0⟩`, as a `Superposition` truncated at `truncation` pairs and trimmed to `threshold`. |
| `transform_state(state, U, V, direction, base_state=None)` | rewrite a `Superposition` in the other basis by repeatedly applying the transformed creation operator to a reference vacuum. |

### State utilities
`normalize_state`, `trim_state`, `get_mode_occupations` (⟨n_i⟩ per mode),
`get_total_excitations`, `partial_trace(state, keep_modes)` (reduced density
matrix, tracing out the complementary modes).

---

## 4. `CompositeEngine`

Bipartite `H_A ⊗ H_B` bookkeeping (here `A` = cavity Fock, `B` = matter Fock):

| method | what it does |
|---|---|
| `vector_to_subsystems(vec, basis_A, basis_B, threshold)` | `1D vector → {state_A: Superposition_B}`. |
| `subsystems_to_vector(state_dict, basis_A, basis_B)` | inverse. |
| `map_subsystem_B(state_dict, map_func)` | apply `map_func` to every `Superposition_B`. |
| `partial_trace(vec, dim_A, dim_B, keep='A'\|'B')` | reduced density matrix of the kept subsystem. |
| `embed_matrix_A(M_A, dim_B)` / `embed_matrix_B(dim_A, M_B)` | `M ⊗ I` / `I ⊗ M`. |
| `kron_matrices(M_A, M_B)` | `M_A ⊗ M_B`. |

---

## 5. Module-level helpers

| function | what it does |
|---|---|
| `compute_position_exponential_sparse(N, u)` | matrix of the displacement operator `exp(i u (a + a†))` in the `N`-level Fock basis (Laguerre formula). Used to build the exact `sin(γ(a+a†))` cavity coupling. |
| `compositions(n, total)` | integer compositions of `total` into `n` non-negative parts (bosonic Fock enumeration). |
| `_eye`, `_destroy`, `_num` | sparse `I`, `a`, `a†a` on a truncated Fock space. |

---

## 6. `PHOTONICChain`

```python
PHOTONICChain(N, Omega_C, Chi_C, omega_r, gamma, fock_photon,
              Omega_J, Chi_J, kerr, PBC=False,
              coeff_tol=1e-9, zero_tol=1e-12, statistics='boson')
```

| arg | meaning |
|---|---|
| `N` | unit cells; `n = 2N` sites |
| `Omega_C`, `Chi_C` | capacitive coupling scale and SSH dimerization |
| `Omega_J`, `Chi_J` | inductive (Josephson) coupling scale and dimerization |
| `omega_r` | cavity frequency (qubit-frequency units) |
| `gamma` | qubit↔cavity coupling (argument of the `sin`) |
| `fock_photon` | cavity Fock-space truncation |
| `kerr` | on-site anharmonicity magnitude (`0`/`None` → no Kerr; forced for `statistics='fermion'`) |
| `PBC` | periodic vs open chain |
| `statistics` | `'boson'` or `'fermion'` (fermions: no Kerr, antisymmetric pairing) |

### What the constructor does

1. Builds the capacitive (`A, B`) and Josephson (`AJ, BJ`) Nambu blocks.
2. **Bogoliubov-diagonalizes the capacitive quadratic form** → quasiparticle
   ("eta" / "new") modes with frequencies `ω_j`, and rotates the Josephson
   sector into that basis.
3. **If `kerr ≠ 0`:** a two-pass procedure — diagonalize, transform the Kerr
   tensor into the eta basis and normal-order it, fold the residual quadratic
   `(1,0)` piece it produces back into the capacitive form, then re-diagonalize.
   `self.kerr_operator` ends up as the `(1,1,0,0)` `TensorTerm` `V_{ijkl}` in the
   *final* eta basis.

### Attributes set

| attribute | meaning |
|---|---|
| `n` | `2N` |
| `U`, `V` | `(n,n)` Bogoliubov coefficients, site → eta |
| `P`, `inverse_P` | `(2n,2n)` pseudo-unitary `[[U,V*],[V,U*]]` and its inverse |
| `H_C` | `(2n,2n)` capacitive Hamiltonian in the eta basis (diagonal). Physical frequencies are `_nambu_prefactor() · diag(H_C)[:n]` (×2 bosons, ×1 fermions). |
| `H_J` | `(2n,2n)` Josephson Hamiltonian in the eta basis (**not** diagonal) |
| `kerr_operator` | `TensorTerm` `V_{ijkl}` `(1,1,0,0)` in the eta basis, or `None` |
| `engine` | the `SecondQuantizationEngine` (`n` modes) |

---

### 6.1 Fixed-excitation Fock-space Hamiltonians

| method | what it does / physical scope |
|---|---|
| `matter_hamiltonians(n_excitations, returns_subspace=False)` | Generate the eta-mode Fock basis with `n_excitations` quanta and build the matter Hamiltonian there: `H_C_fock` (diagonal eta energies; if Kerr and `n_excitations ≥ 2`, adds the Kerr matrix and diagonalizes it, storing the within-sector rotation as `self.transform`) and `H_J_fock` (number-conserving Josephson part, rotated by the same `transform`). With `returns_subspace=True` returns everything instead of storing it. Sets `fock_states`, `fock_state_index`, `fock_dim`, `fock_energies`, `H_C_fock`, `H_J_fock`. |
| `build_HC_fock(states, n_total, save=True)` | Diagonal capacitive Hamiltonian on an explicit Fock basis (the Kerr-diagonalization branch of `matter_hamiltonians`). |
| `build_HJ_fock(states, state_index, n_total, transform=None)` | Number-conserving Josephson Hamiltonian `Σ (A_J + D_J^T)_{ij} b†_i b_j` on an explicit Fock basis, optionally conjugated by `transform`. |
| `set_basis(states)` | Register an explicit matter Fock basis (`fock_states`, `fock_dim`, `fock_state_index`) without building any Hamiltonian. |
| `single_particle_translation_op()` | Build the single-particle lattice translation operator and its eta-basis, normal-ordered form (`self.translation_op`, `self.transformed_translation_op`) — for crystal-momentum / band analysis. |

### 6.2 Full cavity + matter Hamiltonian and spectrum

| method | what it does |
|---|---|
| `build_full_hamiltonian(linear_coupling=True, sw=False)` | Assemble `self.H` (sparse) on `photon ⊗ matter` for the sector set by `matter_hamiltonians`. `linear_coupling=False` uses the exact `sin(γ(a+a†))`. `sw=True` adds the Schrieffer–Wolff correction. Also stores `self.sin_op`. |
| `schrieffer_wolff_correction(eps_denom=1e-6)` | Second-order SW that perturbatively eliminates the ±2-excitation (pair) Josephson terms coupling matter to the cavity, giving an effective in-sector correction `self.H_SW`. Requires `sin_op`, `fock_energies`. |
| `get_low_lying_spectrum(k=5)` | lowest `k` eigenpairs of `self.H` (`scipy.sparse.linalg.eigsh`, Hermitized). |

### 6.3 Probe-resonator dispersive shift

**Scope.** Given the *already-diagonalized* chain (eta modes `ω_j`, Kerr tensor
`V`, assumed decoupled from the original cavity), attach a **new probe
resonator** `α` of frequency `omega_r_new` coupled to the sites as
`(α† − α) Σ_j g_j (b†_j − b_j)` and compute the **state-dependent dispersive
shift** `χ_i` the probe acquires for each single excitation of the chain. Only
meaningful in the dispersive regime; the shift is driven entirely by the Kerr
nonlinearity, so it is identically zero when `kerr = 0`.

**Pipeline.** Transform `g_j → g'_j` (eta basis) → Bogoliubov-diagonalize the
`(n+1)`-mode quadratic part into modes `β_k` → pick `β_0`, the one closest in
frequency to `omega_r_new` → push `V` through that transform → keep only the
`β_0†β_0`-weighted quadratic chain piece `(Σ_ij V''_ij β†_i β_j) β_0†β_0` → block-
diagonalize `V''` on degenerate chain subspaces → second order:
`χ_i = V''_ii + Σ_{j≠i} |V''_ij|² / (E_i − E_j)`.

| method | what it does |
|---|---|
| `compute_coupling_vec_from_sites(g_sites)` | transform an arbitrary length-`n` per-site coupling vector into the eta basis → `g'` (length `n`). |
| `compute_coupling_vec(g0, ground_mode, mode, positions)` | same, for the standing-wave parametrization `g_sites[j] = g0 · sin(2π(1+mode) x_j / (L·ground_mode))`. `positions` are normalized `x_j / L`. |
| `compute_resonator_mode_transformation(omega_r_new, g0, ground_mode, mode, positions)` | run only the second Bogoliubov diagonalization; returns `(omega_beta, U2, V2, P2, b0, g_prime)`. |
| `transform_kerr_to_resonator_basis(U2, V2)` | push `self.kerr_operator` through the `(n+1)`-mode transform, normal-order, return the `(1,1,0,0)` tensor `V'` (shape `(n+1)⁴`). Raises if no Kerr. |
| `compute_resonator_dispersive_shift_from_gprime(omega_r_new, g_prime, degeneracy_tol=1e-6, dispersive_tol=0.1)` | **the shared core** — the whole pipeline, keyed on the already-transformed `g'`. Returns the length-`n` real vector `χ`, indexed by ascending `β` frequency (probe removed), after the degenerate-subspace rotation. Raises `ValueError` if any `|g'_j| / |ω_j − omega_r_new| ≥ dispersive_tol`. |
| `compute_resonator_dispersive_shift_vec(omega_r_new, g0, ground_mode, mode, positions, …)` | standing-wave entry point → core. |
| `compute_resonator_dispersive_shift_vec_from_sites(omega_r_new, g_sites, …)` | arbitrary per-site entry point → core. |

*Accuracy:* dropping every term outside the chain's single-excitation subspace
makes the result exact to `O(kerr)` — residual `O(kerr²)`, i.e. a few percent
when `kerr` is comparable to the chain level spacing, vanishing as `kerr`
shrinks (cross-checked against a full QuTiP diagonalization).

Private helpers: `_nambu_prefactor`, `_extended_engine`,
`_diagonalize_extended_quadratic`, `_resonator_mode_transformation_from_gprime`,
`_resonator_number_kerr_matrix` (extracts `V''`),
`_diagonalize_degenerate_blocks`.

**Mode-discriminability tooling** — how well a single dispersive measurement
separates the chain modes:

| name | what it does |
|---|---|
| `probe_freq_in_gap(gap_index, fraction=0.5, margin=0.4)` (method) | an `omega_r_new` placed relative to the chain's sorted eta spectrum: `gap_index = -1` below every mode, `k ∈ 0..n-2` inside the gap between sorted modes `k, k+1`, `n-1` above every mode. Interior gaps land at `fraction` through the gap; the below/above positions sit `margin × bandwidth` past the extreme mode (a bandwidth fraction, not a spacing count — a uniform `g_j` couples mostly to one band edge, so "above all" would otherwise be near-resonant with it). Because the spectrum moves with `(Omega_C, Chi_C)`, calling this per grid point keeps "the probe is in gap k" true across a sweep. |
| `mode_discriminability(omega_r_new, g_sites, normalize=1.0, …)` (method) | `min_{i<j}|χ_i − χ_j| / normalize` — the smallest pairwise gap between the per-mode shifts (can I tell two excited modes apart). Pass `normalize = g**2` for a coupling-independent figure of merit (`χ ~ O(g² kerr)`). Propagates the dispersive-regime `ValueError`. |
| `min_abs_shift(omega_r_new, g_sites, normalize=1.0, …)` (method) | `min_i |χ_i| / normalize` — the weakest per-mode shift (does every mode leave a visible imprint at all). Same normalization/error behaviour. On a grid it also drops straight out of `sweep_discriminability`'s `chi_cube`: `np.min(np.abs(chi_cube), axis=-1)`. |
| `min_pairwise_gap(values)` (module function) | `min(diff(sort(values)))`; `+inf` for < 2 values. |
| `sweep_discriminability(make_chain, omega_C_values, chi_C_values, g_sites, *, gap_index=None, omega_r_new=None, normalize=1.0, …)` (module function) | build `make_chain(Omega_C, Chi_C)` at every grid point, place the probe (adaptive `gap_index`, or fixed `omega_r_new`), record `min_pairwise_gap(χ)/normalize`; non-dispersive points → `NaN`. Returns a dict: `gap_grid`, `chi_cube` (…, 2N), `omega_r_new_grid`, `omega_C_values`, `chi_C_values`. |

### 6.4 Single-particle dispersive shift and common bus (linear, no Kerr)

| method | what it does / scope |
|---|---|
| `compute_dispersive_shift_vec(Omega, g0, ground_mode, mode, positions, dispersive_tol=0.1)` | second-order shift `g'_j² / (ω_j − Ω)` between one resonator mode (frequency `Omega`, standing-wave profile) and each chain eta mode. Raises if `|g'_j|/|ω_j − Ω| ≥ dispersive_tol`. Flip the overall sign for the resonator-minus-qubit convention. |
| `compute_common_bus_shift_matrix(Omega, positions, g0=1, ground_mode=2, dispersive_tol=0.1)` | a whole multimode bus: bus mode `m` at frequency `Omega·(m+1)` (harmonic ladder), capacitive coupling `g0·√(m+1)`. Returns the `(n × n)` matrix of dispersive shifts (bus mode × chain mode). |
| `compute_common_bus_coupling_matrix(positions, g0=1, ground_mode=2)` | same bus, but the `(n × n)` matrix of bare eta-basis couplings `g'` (same `√(m+1)` scaling). |

### 6.5 Vacua and basis transforms (bare ↔ Bogoliubov)

The Bogoliubov ("new"/eta) vacuum is a squeezed state in the bare ("old") Fock
basis, and vice versa.

| method | what it does |
|---|---|
| `build_squeezed_vacuum(truncation, threshold=0.99)` | eta-vacuum written in the bare basis → stores `self.vacuum` (dict). |
| `build_unsqueezed_vacuum(truncation, threshold=0.99)` | bare vacuum written in the eta basis → stores `self.unsqueezed_vacuum`. |
| `new_to_old_fock(new_fock_state, threshold=0.99, use_cache=True)` | expand one eta-basis Fock state in the bare basis (needs `build_squeezed_vacuum` first). Cached. |
| `old_to_new_fock(old_fock_state, threshold=0.99, use_cache=True)` | expand one bare Fock state in the eta basis (needs `build_unsqueezed_vacuum` first). Cached. |
| `build_new_to_old_matrix(truncation=8, threshold=0.99)` | assemble the (rectangular) basis-change matrix for the current `fock_states` sector, new → old. Returns `(matrix, row_states, row_index, col_states)`. |
| `build_old_to_new_matrix(old_states, …)` | same, old → new, for an explicitly supplied list of old states. |
| `build_squeezed_vacuum_tensor(alpha=0, …)` / `build_unsqueezed_vacuum_tensor(alpha=0, …)` | the respective vacuum as a full `photon ⊗ matter` composite vector, with the cavity in Fock state `alpha`. |
| `transform_tensor_to_old(full_vector, truncation=2, threshold=0.99)` / `transform_tensor_to_new(full_vector, truncation=8, threshold=0.99)` | transform an arbitrary composite state vector between bases. |
| `warm_old_basis_cache(truncation=8, threshold=0.99)` | precompute the new→old matrix (`self._C_matrix`) used by `compute_ipr_batch`. |

### 6.6 State decomposition, partial trace, observables

| method | what it does |
|---|---|
| `decompose_tensor_state(full_vector, threshold=0.99)` | composite vector → `{photon_number: {matter_fock_state: amplitude}}`, each matter block trimmed to `threshold`. |
| `compose_tensor_state(decomposition)` | inverse; returns `(psi, basis_B, fock_state_index)`. |
| `decompose_to_old_fock(full_vector, threshold=0.99)` | as `decompose_tensor_state` but with each matter block further expanded into the bare Fock basis. |
| `partial_trace(full_vector, which_to_keep='matter')` | reduced density matrix of `'matter'` or `'photon'`. |
| `quanta_distribution(fock_vector, …)` | normalized per-site excitation profile `⟨n_j⟩` of a matter state (accepts a vector over `fock_states` or a `{state: amp}` dict). |
| `quanta_ipr(fock_vector)` | inverse participation ratio of that site-occupation profile (localization measure); returns `(ipr, 1/ipr, distribution)`. |
| `quanta_ipr_mixed(rho_matter)` | same from a matter density matrix. |
| `compute_ipr_batch(vecs)` | batched IPR over many full eigenvectors at once (requires `warm_old_basis_cache` first). |

---

## 7. Typical workflows

**Band structure / quasiparticle spectrum of the bare chain**
```python
chain = PHOTONICChain(N=6, Omega_C=0.6, Chi_C=0.3, omega_r=0.1, gamma=0.0,
                      fock_photon=1, Omega_J=1.0, Chi_J=0.0, kerr=0.0)
omega_eta = 2.0 * np.diag(chain.H_C)[:chain.n].real     # eta-mode frequencies
```

**Full cavity + chain Hamiltonian in the 1-excitation sector**
```python
chain.matter_hamiltonians(n_excitations=1)
chain.build_full_hamiltonian(linear_coupling=True, sw=False)
evals, evecs = chain.get_low_lying_spectrum(k=10)
```

**Dispersive readout of a probe resonator (needs Kerr)**
```python
chain = PHOTONICChain(N=3, Omega_C=0.6, Chi_C=0.35, omega_r=0.1, gamma=0.01,
                      fock_photon=6, Omega_J=1.0, Chi_J=0.0, kerr=0.03)
chi = chain.compute_resonator_dispersive_shift_vec_from_sites(
    omega_r_new=0.62, g_sites=np.full(chain.n, 3e-3), dispersive_tol=0.3)
```
