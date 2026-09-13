# Flux-coupling constant `c` for the SQUID-array design — extraction log

**Design:** `build_design.py` — a chain of `N` series SQUID loops (`SquidArray`,
no transmon) running down the centreline of a hairpin fold of a meandered CPW
readout resonator (`readout_1`), on the CNEA cluster
(`arturo.ruiz@10.73.25.223`, SGE).

**Bottom line up front:** `c ≈ 2.3e-3 rad` for `n = 23`, `L_J = 600 nH`, from
the direct B-field method (§3) — the flux coupling that modulates each SQUID
phase. Per-loop mutual `M_k ≈ 1 pH` (0.98–1.72 pH along the standing wave),
**flat**; collective `Σ M_k ≈ 21–27 pH`; `L_r ≈ 0.87 nH`. The coupling is
**distributed** — every loop contributes `c_k ≈ 9e-5 rad`, same sign,
`c_array = Σ c_k` **linear in n**. The port-EPR total (§4) is ~4.7e-3, higher
because the chain-terminal junctions carry ~10× the interior current (current
distribution, not flux). Bracket `c ≈ 2–5e-3`. Shape robust; prefactor to ~2×.

**The method is normalisation-independent.** A Palace eigenmode has an
arbitrary overall scale `s`; `c`, `M`, `L_r` are ratios (`Φ/√U_m`, `Φ/I_r`) and
are `s`-invariant — verified: rescaling B,E by `s` and U_m by `s²` for
`s = 1e-6…1e3` moves them by `<1e-6`. Palace's field *output* is SI (docs:
"all fields written out in SI units; post-processing mesh in units of L0 m");
confirmed 3 ways (`∫U_m dV = domain-E E_mag` to 5 figs; `U_m` density
`= (|Br|²+|Bi|²)/(2μ0)` exactly; port `|V/I| = |Z(L_J‖C_J)|`).

**Two bugs, both fixed 2026-08-29:**
1. *Mesh length unit.* `L0 = 1e-3` → mesh coords in **mm**; pyvista areas in
   mm², volumes in mm³ while B, U_m are SI. Fix: `area×L0²`, `vol×L0³`. Without
   it `Φ×1e6`, `U_m_tot×1e9`, `M×1e6`, `c×31.6` (this was the `c≈0.05, M≈1µH`
   error).
2. *Energy convention.* Palace's `U_m` density is `(|Br|²+|Bi|²)/(2μ0)` = the
   **peak** inductive energy `½L_r I_r²`, not the time average. Fix:
   `L_r = 2 U_m_tot/I_r²` (not 4), `c_k` denominator `4 U_m_tot` (not 8) — a
   `√2` in `c`. `M_k` unaffected. `L_r ≈ 0.87 nH` now matches a hand estimate.

---

## 1. Why L_J was raised to 600 nH

At the original `L_J = 15 nH` the SQUID-array "phonon" band (the array is a
~1.1 mm slow-wave line) spans ~3–13+ GHz and straddles the ~8 GHz readout
resonator — the resonator mode sat 0.25 GHz from an array phonon and
hybridised. Band top ~ `1/sqrt(L_J)`:

| L_J | phonon band top | gap to ~8 GHz resonator |
|---|---|---|
| 15 nH | ~30 GHz | overlaps |
| 400 nH | ~5.8 GHz | ~2 GHz |
| 600 nH | ~4.7 GHz | ~3 GHz — clean |

Confirmed by a band-check run: at L_J = 600 every array mode solved from
3.5 GHz came back below 5.05 GHz. Raising `L_J` leaves the CPW resonator
frequency essentially unchanged (`f_r ≈ 7.94 GHz`, L_J-independent to 4 dp).
The dispersion capacitance is the loop-metal-to-ground `C_g` (3 µm gap), not
the 1.5 fF junction `C_J`.

---

## 2. Mesh / AMR

AMR hit the memory wall: 46 lumped ports + AMR's retained meshes / error
estimation / rebalancing is ~3× a plain solve; the (600, 23) convergence AMR
(`1938535`) reached **501 / 503 GB** on the highmem node and was `qdel`'d
before it could OOM. Frequency was converging (Δf 1.0% → 0.56%) but the port
participation was not — which is now moot: the B-field method does not use the
ports. All results below are **non-AMR fixed-mesh** solves (memory-safe at
~40–150 GB).

Meshes used:
- `fullarray/a1.5/` — whole-array `fine_mesh_components` at 1.5 µm min,
  **3.77 M elements**, `f_r = 7.9355 GHz`. The headline mesh.
- `final/n{5,10,16,23}_box2.0/` — coarser (`array-bg-min 4 µm`), `f_r ≈ 7.946`.
  The length sweep.

---

## 3. The B-field method (the primary result) — `extract_M_bfield.py`

Derivation and assumptions: `C_CALCULATION.md` §2b and the script docstring.

```
c_k = π Φ_k sqrt(ħ ω_r / (4 U_m_tot)) / Φ0        (needs only Φ_k and U_m_tot — no current model)
M_k = Φ_k / I_r                                    (needs the CPW-current model; ~4× uncertain)
L_r = 2 U_m_tot / I_r^2                            (U_m_tot = peak E_mag, per domain-E.csv)
```

`Φ_k` = ∫ B_z over loop k's interior on a cut a few µm above the metal plane
(B_z kept complex — eigenmode phase is in B_imag, `|B_imag|/|B_real| ~ 8e5`);
`U_m_tot` = ∫ U_m dV = Palace's `E_mag`.

### 3.1 Self-consistency

- `M_k^loc = |Φ_k| / I_local` (flux over the *local* resonator current) is
  **flat to ~11%** across all loops, all n, all flux-cut heights — geometric.
- `c_k` direct and via `π M_k I_r^zpf / Φ0` are **identical**.
- **normalisation-invariance:** rescale B,E by `s`, U_m by `s²`, `s = 1e-6…1e3`
  → `c, M, L_r` move `<1e-6`.
- `∫ U_m dV_SI` = `domain-E.csv E_mag = 5.90408e-12 J` to 5 figures;
  `U_m` density `= (|Br|²+|Bi|²)/(2μ0)` exactly; port `|V/I| = 2.418e4` =
  `|Z(L_J‖C_J)| = 2.42e4`.

### 3.2 Per-loop profile (n = 23, L_J = 600, a1.5 mesh, z = 3 µm)

`Φ_k` a smooth ramp, loop 0 `2.0e-13` → loop 22 `1.14e-13` Wb (factor ~1.8).
`M_k` runs `1.72 → 0.98 pH` — **flat**. `c_k` runs `1.44e-4 → 8.2e-5`, all
same sign. `L_r = 0.87 nH`.

```
c_array (signed = Σ|c_k|) = 2.28e-3 rad
Σ M_k = 27.2 pH   (collective array-to-resonator mutual)
```

### 3.3 c_array and Σ M_k vs number of loops (z = 2.5 µm, box2.0 mesh)

| n | f_r (GHz) | c_array (rad) | Σ M_k (pH) | c_array / n |
|---|---|---|---|---|
| 5  | 7.9461 | 4.52e-4 | 5.0  | 9.0e-5 |
| 10 | 7.9402 | 8.74e-4 | 9.6  | 8.7e-5 |
| 16 | 7.9453 | 1.38e-3 | 15.8 | 8.6e-5 |
| 23 | 7.9461 | 1.92e-3 | 21.1 | 8.4e-5 |
| 23 (a1.5 mesh) | 7.9355 | 2.28e-3 | 27.2 | 9.9e-5 |

**Linear**, slope `~9e-5 rad` per loop. a1.5 gives ~19% more than box2.0 —
finer = larger, not fully plateaued.

### 3.4 c_array vs flux-cut height (a1.5, n = 23)

| z (µm) | c_array (rad) | Σ M_k (pH) |
|---|---|---|
| 2.0 | 2.221e-3 | 25.8 |
| 2.5 | 2.249e-3 | 26.6 |
| 3.0 | 2.276e-3 | 27.2 |

Stable to ~3% over 2–3 µm. Quote **`c ≈ 2.3e-3 rad`** (a1.5) / `1.9e-3` (box).

### 3.5 The junction-corner field spike

`render_bfield_maps.py`. The user's ParaView close-up showed `B_z` 2–3 orders
above the surrounding loop field, in many loops. It is a **sharp-corner (PEC)
singularity** at the 90° edges where the junction-port surfaces meet the loop
plates (`y ≈ ±8–15 µm`, on the loop-cell boundaries). Bright "bowtie" features
at `z = 0.3–1 µm`; smoothed away by `z = 2.5 µm`.

**It does not affect the flux integral.** The window (`|y| < 12.5 µm`, inset
from the plates) excludes the corners. Dropping the brightest 1% of cells
inside the window changes `Φ_k` by **< 2%** at every height; max/median in the
window is 1.2–2.4. Checked for loops 3, 11, 20.

### 3.6 Meissner / PEC

Palace enforces `B·n̂ = 0` on every surface and carries the surface current
`J_s` — that **is** the Meissner effect at the level of perfect screening. The
array body + the 3 µm-away ground plane screen the resonator flux, and the
result includes it (`M_k^loc ≈ 1 pH`, below the bare geometric value).

**Not modelled:** the London depth `λ_L ~ 50–100 nm` / kinetic inductance
(SQDMetal has no `λ_L`). It folds into an effective `L_loop` added to
`2 L_J = 1200 nH`; at `L_kin ~ pH` a `< 1e-5` correction to `c`. The `L_J`
= 400-vs-600 check (§4.1) bounds it. DC Meissner / trapped flux does not
enter — `c` is the small-signal fluctuation about the DC bias and the integer
in `Φ = nΦ0` differentiates away.

### 3.7 Raw data

`bfield_data/*.npz` — one file per (n, mesh, z). Keys (all SI): `loop`,
`Phi_complex`, `Phi_signed`, `Phi_abs`, `I_cpw`, `I_twin`, `I_r_peak`,
`I_r_zpf`, `M_k`, `M_loc`, `c_k`, `c_k_viaM`, `p_EPR`, `p_pred_from_M`,
`L_r_H`, `U_m_tot_J`, `E_mag_domainE_J`, `f_r_Hz`. Load with `numpy.load()`.

---

## 4. The port-EPR path (independent cross-check) — `extract_c_array.py`

Reformulation of `../squid_resonator_epr` for an array: each loop's two
junctions both carry `L_J`, so the loop is a `2 L_J` circuit,
`p_loop_k = p_{2k} + p_{2k+1}`, and

```
c_k = (1/2) sqrt( |p_loop_k| f_r L_J[nH] / 163.4 )
```

The port-EPR sum is dominated by loops 0, 21, 22 and gives `c ≈ 3–5e-3`.

### 4.1 Three Palace outputs per loop (`audit_units_and_currents.py`, job 1938551)

| quantity | interior loops | end loops (0, 21, 22) |
|---|---|---|
| `M_k` from B-field flux | ~1.1 pH | 0.98–1.72 pH |
| `M_k` from `sqrt(2 L_J L_r p_EPR)` | ~1.8 pH (agree ~1.5×) | 11–46 pH (disagree ~30×) |
| junction current (`port-I.csv`) | ~1.5e-7 A | 1–4e-6 A |
| within-loop junction-current asymmetry | 0.3–0.5 | 0.03–0.13 |

**The interior loops now agree** between the two methods (`M_k ~ 1 pH` both
ways, ~1.5×). The **end loops** disagree ~30×. `port-I.csv` explains it: the
chain-terminal junctions carry ~10× the interior current, and their two
junctions carry *equal* current (asymmetry → 0) while interior loops show a
30–50% within-loop mismatch. The resonator-induced current **circulates
cleanly through the terminal loops** (the chain ends there, no other path) but
**partly bypasses the interior junctions** via the link bars / continuous
array body. That extra terminal current is inductive-energy participation, not
loop-threading flux — it lifts the port-EPR total to ~4.7e-3 but does not
change the geometric `M_k`.

The earlier "only the two ends couple — the array screens its interior"
reading is **wrong**: `M_k` (the flux) is flat across all loops. port-EPR's
end weighting is a current-distribution effect at the chain terminals. The
user's ParaView read — "the eigenmode doesn't particularly have more flux on
the ends" — is correct.

**Which `c`:** the flux-coupling `M_k` (→ `c ≈ 2.3e-3`) modulates each SQUID
phase in the Hamiltonian; the port-EPR total (~4.7e-3) is closer to the
array's inductive loading of the resonator. Bracket `c ≈ 2–5e-3`.

### 4.2 L_J two-point check (400 vs 600 nH)

Via port-EPR on the end loops, `c(400)` and `c(600)` agree to ~5–6% at every
mesh → `L_loop << 2 L_J`, so `c ≈ c_∞`, no screening-inductance fit. This is
the experimental bound on the un-modelled kinetic inductance (§3.6).

---

## 5. Validation checklist (EPR_C_EXTRACTION.md §6)

- [x] **6.1 two-point L_J** — `c(400) ≈ c(600)` to ~6% → `c ≈ c_∞`, `L_loop`
  and `L_kin` negligible.
- [~] **6.2 junction-field / c convergence** — the B-field method sidesteps the
  port-participation floor. `c_array` stable to ~3% over flux-cut heights
  2–3 µm; ~18% mesh sensitivity (box vs a1.5), finer = larger, not fully
  plateaued.
- [x] **6.3 mode ID by participation, not frequency** — `identify_resonator_mode`,
  resonator = lowest-f mode with `Σ|p_j| << 1`.
- [x] **c < 0.1 sanity** — per loop `c_k ≈ 9e-5`; `c_array ≈ 2.3e-3` coherent
  over 23 loops.
- [x] **units** — SI (Palace docs + 3 checks); `∫ U_m dV = domain-E E_mag` to
  5 figs; `U_m` density `= (|Br|²+|Bi|²)/(2μ0)`; port `|V/I| = |Z(L_J‖C_J)|`.
- [x] **normalisation-invariance** — rescale fields by `s = 1e-6…1e3` → `c, M,
  L_r` move `<1e-6`.
- [x] **method self-check** — `M_k^loc` flat to ~11%; `c_k` direct == `c_k` via M.
- [~] **B-field vs port-EPR** — agree on the interior (~1 pH both); disagree
  ~30× on the end loops (chain-terminal current distribution, §4.1). `c`
  bracket: 2–5e-3.

---

## 6. Files

| file | role |
|---|---|
| `extract_M_bfield.py` | **the B-field method** — Φ_k, CPW-current inversion, M_k, c_k, L_r; writes the .npz |
| `audit_units_and_currents.py` | units / normalisation-invariance / port-I audit (§3.1, §4.1) |
| `bfield_data/*.npz` | raw per-loop arrays (all SI); one per (n, mesh, z) |
| `render_bfield_maps.py` | |B| / B_z colormaps over a loop + singularity-vs-height check |
| `build_design.py` / `build_and_eigenmode.py` | geometry + eigen solve, 2N lumped ports |
| `extract_c_array.py` | port-EPR path — the independent cross-check |
| `C_CALCULATION.md` | full derivation, assumptions, validation |

Cluster runs under `~/circuit_design/squid_array_coupling/`:
`fullarray/a1.5/`, `final/n{5,10,16,23}_box2.0/`, `ljscan/`. B-field re-run
(corrected L_r): job 1938552; units audit: job 1938551; `.npz` in each run dir
and in `bfield_data/`.
