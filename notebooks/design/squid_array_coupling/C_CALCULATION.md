# The `c` calculation for the SQUID array — full spec

This is the `c` calculation, written out so you can check it before trusting
any number. Two paths: the **direct B-field method** (`extract_M_bfield.py`,
§2b) is the one used for the reported result; the **port-EPR path**
(`extract_c_array.py`, §3) is a cross-check that only resolves the array ends
(§7.1). It is a reformulation of
`../squid_resonator_epr/EPR_C_EXTRACTION.md` and
`../squid_resonator_epr/C_PARAMETER_DERIVATION.md` for an array of loops
instead of one loop. Read those two first; this doc only covers what changes.

Constants: `phi0 = Phi0/2pi = 3.29106e-16 Wb`, `h`, and the convenient combo
`phi0^2 / h = 163.4 GHz * nH`.

---

## 1. What `c` is (unchanged from the single-SQUID case)

The resonator's current threads flux through a SQUID loop. The loop's two
junction phases enter the Hamiltonian as

```
-E_J cos(phi_a) - E_J cos(phi_b) = -2 E_J cos((phi_a+phi_b)/2) cos((phi_a-phi_b)/2)
```

Fluxoid quantization ties the differential phase to the loop flux,
`phi_a - phi_b = 2 pi Phi_loop / Phi0`. Split that flux into a DC bias plus the
resonator's zero-point fluctuation, `Phi_loop = Phi_dc + M * I_r^zpf (a + a^dag)`:

```
(phi_a - phi_b)/2  =  phi_dc/2  +  c * (a + a^dag)
c  =  pi M I_r^zpf / Phi0  =  M I_r^zpf / (2 phi0)
```

`c` is a **linear-circuit quantity** — a mutual inductance `M` (resonator mode
<-> loop) times the resonator's zero-point current, over the flux quantum. It
is the only piece needing a field solver. Expected size: `1e-3` to `1e-2` rad;
`> 0.1` means something is wrong.

---

## 2. What Palace gives you: the participation ratio

You cannot read `M` off an open loop — no closed path, no current, zero
coupling. The Minev-2021 trick: put a **finite lumped inductor at the junction
gap** and let Palace report the fraction of each eigenmode's inductive energy
stored in that lumped port:

```
p  =  (inductive energy in the lumped port)  /  (total inductive energy of the mode)
```

written to `port-EPR.csv`, one column per inductive lumped port, one row per
mode.

**Single SQUID (reference):** one arm is a lumped `L_probe`, the other arm is
shorted with metal, so the loop is a one-inductor circuit and

```
c_m  =  (1/2) sqrt( p_probe * f_m / (2 * E_probe/h) ),   E_probe/h = phi0^2/(L_probe h) = 163.4 / L_probe[nH]  GHz
```

Derivation (see C_PARAMETER_DERIVATION.md 2.1): with `I_loop = M I_r / L_probe`
(the loop flux is absorbed by the probe inductor), `p_probe = M^2/(L_probe L_r)`
where `L_r` is the mode's lumped-equivalent inductance; substitute into
`c^2 = pi^2 M^2 (h f_r / 2 L_r) / Phi0^2` and the `L_r` and `M^2` cancel out,
leaving the formula above.

---

## 2b. The direct B-field method (the one actually used) — `extract_M_bfield.py`

The participation path above **fails for this geometry**: the whole
array<->resonator coupling is a `~1e-6` fraction of the mode energy, and the
per-loop participation of the ~20 interior loops sits at `|p_loop| ~ 1.5e-9`
with random sign — below what the eigensolve resolves against a mode living
almost entirely in the CPW. `port-EPR.csv` only lifts off its numerical floor
for the 3 loops at the ends of the chain. See 7.1.

So `M_k` is read straight off the eigenmode fields instead, with no lumped
port:

1. **Loop flux.** `Phi_k = integral of B_z over loop k's interior`, on a
   horizontal cut a few um ABOVE the metal plane (on `z=0` itself `B.n = 0` by
   the PEC BC). The eigenmode's phase is in the imaginary part here
   (`|B_imag|/|B_real| ~ 8e5`), so `B_z` is kept complex, globally
   phase-aligned (`phi0g = angle(Phi[argmax|Phi|])`), then `Phi_s = Re(Phi
   e^{-i phi0g})`.

2. **Resonator modal current**, from its OWN CPW gap field via the exact
   zero-thickness-PEC result (Lahl & Woerdenweber, APL 81 505 (2002);
   Bothner arXiv:1202.5963, lambda_eff -> 0):

   ```
   center |x|<a, gaps a<|x|<b, k = a/b, K = complete elliptic integral 1st kind
   B_z,gap(x) = mu0 I b / ( 4 K(k) sqrt( (x^2-a^2)(b^2-x^2) ) ),   a<|x|<b
   ```

   inverted at a sampled gap point. Our resonator: `a = 10 um`, `b = 22.25 um`,
   `k = 0.449`, `K(k) = 1.66`.

3. **Normalisation.** `L_r = 2 U_m_tot / I_r^2`, `I_r^zpf^2 = hbar w_r /
   (2 L_r)`, and `I_r` cancels:

   ```
   c_k = pi * Phi_s,k * sqrt( hbar w_r / (4 U_m_tot) ) / Phi0
   ```

   so **`c_k` needs only `Phi_k` and `U_m_tot`** -- NOT the current model.
   `M_k = Phi_s,k / I_r` does, and inherits the current model's ~4x
   uncertainty (exact isolated-CPW vs twin-lead give `I_r` differing ~4x).
   `L_r ~ 0.87 nH` (matches a hand estimate for a lambda/2 CPW antinode,
   L' l / 2).

4. **Total.** `c_array = sum_k c_k`, `M_array = sum_k M_k` (all loops thread
   the resonator flux the same way -> same sign -> coherent sum). No noise
   filter: `Phi_k` is a smooth monotonic ramp, every loop well above zero.

**NORMALISATION INDEPENDENCE.** A Palace eigenmode is defined only up to an
overall scalar `s` (the solver fixes unit coefficient-vector norm). Everything
here is a ratio: `Phi_k ~ s`, `I_r ~ s`, `U_m_tot ~ s^2`, so
`M_k = Phi_k/I_r`, `L_r = 2 U_m_tot/I_r^2` and
`c_k = pi Phi_k sqrt(hbar w/4 U_m_tot)/Phi0` are all `s^0`. Verified: rescaling
B,E by `s` and U_m by `s^2` for `s = 1e-6 .. 1e3` moves `c, M, L_r` by
`< 1e-6`. The absolute field values are NOT used.

**UNITS -- two things that DO matter (both bugs, fixed 2026-08-29):**

- *Mesh length unit.* Palace `config["Model"]["L0"] = 1e-3` -> mesh coordinates
  in MILLIMETRES, while `B` (T) and `U_m` (J/m^3) are SI (Palace docs: "all
  fields are written out in SI units; the post-processing mesh has units of
  L0 m"). pyvista `compute_cell_sizes` returns area in mm^2, volume in mm^3.
  Must apply `area_SI = area_mesh * L0^2`, `vol_SI = vol_mesh * L0^3`. Without
  it `Phi` x1e6, `U_m_tot` x1e9, `M` x1e6, `c` x31.6.
- *Energy convention.* Palace's `U_m` density is `(|Br|^2+|Bi|^2)/(2 mu0)`
  (verified: `[(|Br|^2+|Bi|^2)/(4 mu0)] / U_m = 0.5` exactly), so
  `U_m_tot = sum(U_m dV) = E_mag` in `domain-E.csv` is the PEAK inductive
  energy `(1/2) L_r I_r^2`, not the time average. Hence `L_r = 2 U_m_tot/I_r^2`
  (not 4), and `c_k` has `/(4 U_m_tot)` (not 8) -- a factor `sqrt(2)` in `c`.

Cross-checks that Palace's SI output is real and self-consistent:
`sum(U_m dV_SI) == domain-E E_mag` to 5 figs; port `|V|/|I| = 2.418e4` matches
`|Z(L_J || C_J)| = 2.42e4` for the `L_J=600 nH, C_J=1.5 fF` I set.

**Result:** `c ~ 2.3e-3 rad` for `n=23`, `L_J=600` (B-field flux coupling, all
loops), bracket `2-5e-3` (port-EPR total ~4.7e-3, dominated by the array-end
junctions -- see 7.1). `c_array ~ n * 9e-5 rad` (linear). Per loop
`c_k ~ 9e-5`, `M_k ~ 1 pH` (0.98-1.72 pH); collective `sum M_k ~ 21-27 pH`;
`L_r ~ 0.87 nH`. See C_ARRAY_RESULTS.md §3-4.

---

## 3. What changes for the array

`SquidArray` carves a real gap in **both** junctions of **every** loop and puts
`L_J` in each — nothing is shorted. So for loop `k`:

### 3.1 The loop is a `2 L_J` circuit, not an `L_probe` circuit

The circulating current `I_k` driven by the resonator flux goes **through both
junctions in series** around the loop:

```
i w M_k I_r  =  i w (L_J + L_J) I_k     ==>     I_k = M_k I_r / (2 L_J)
```

(geometric loop self-inductance `L_loop ~ tens of pH` is dropped here; that
approximation is exactly what the L_J sweep in 5 tests.)

### 3.2 The loop's participation is the sum of its two junction ports

```
p_loop,k  =  p_{2k} + p_{2k+1}        (SquidArray emits junction 2k = loop k upper, 2k+1 = lower)
```

and the energy those two ports hold is

```
U_k  =  (1/2) L_J I_k^2  +  (1/2) L_J I_k^2  =  L_J I_k^2
```

so, with `U_mode = (1/2) L_r I_r^2`,

```
p_loop,k  =  U_k / U_mode  =  L_J I_k^2 / U_mode
          =  L_J (M_k I_r / 2 L_J)^2 / ((1/2) L_r I_r^2)
          =  M_k^2 / (2 L_J L_r)                                    ... (*)
```

### 3.3 Per-loop `c_k`

`c_k = pi M_k I_r^zpf / Phi0`, `I_r^zpf^2 = h f_r / (2 L_r)`, `Phi0^2 = 4 pi^2 phi0^2`:

```
c_k^2  =  pi^2 M_k^2 (h f_r / 2 L_r) / (4 pi^2 phi0^2)
       =  M_k^2 h f_r / (8 L_r phi0^2)
```

substitute `M_k^2 = 2 L_J L_r p_loop,k` from (*):

```
c_k^2  =  2 L_J L_r p_loop,k * h f_r / (8 L_r phi0^2)
       =  L_J p_loop,k h f_r / (4 phi0^2)
       =  p_loop,k f_r / ( 4 * phi0^2/(L_J h) )
       =  p_loop,k f_r / ( 4 * 163.4/L_J[nH] )          [ f_r in GHz ]
```

so

```
c_k  =  (1/2) sqrt( p_loop,k * f_r * L_J[nH] / 163.4 )
     =  (1/2) sqrt( p_loop,k * f_r / (2 * E_eff/h) ),    E_eff/h = 163.4 / (2 L_J[nH])  GHz
```

**This is exactly the reference formula with `L_probe -> 2 L_J` and
`p_probe -> p_loop,k`.** The `/4` vs the reference's `/8` and the `2 L_J` vs
`L_probe` both come from (*)'s factor of 2, which is the "two junctions, not
one arm + a short" difference.

The sign of `c_k` is carried from `sign(p_loop,k)` — Palace's per-port
participation is signed, and the sign tracks `sign(M_k)`, i.e. which way the
resonator's modal current threads that loop.

### 3.4 Total array coupling

Each loop contributes a term `c_k (a + a^dag)` to its own differential phase.
For the array acting as one collective flux coupler the per-loop mutual
inductances add in series (the resonator flux threads every loop), so
`M_total = sum_k M_k` and, since `c_k` is linear in `M_k`,

```
c_array  =  sum_k c_k                                   (signed sum)
```

Also reported:
- `sum_k |c_k|`  — incoherent upper bound (relevant if the signs are unreliable)
- the dominant single loop's `c_k`

**Noise filter.** Loops that carry no real coupling sit at `|p_loop,k| ~ 1e-9`
with random sign; summed over ~20 loops they can swamp a true `c ~ 1e-3`. So
`extract_c_array.py` keeps only loops with `|p_loop,k| >= 0.08 * max_k |p_loop,k|`
and sums `c_k` over those.

---

## 4. Assumptions (all inherited from the single-SQUID method, plus one)

| # | assumption | why it's OK here |
|---|---|---|
| 1 | loop inductance = `2 L_J`, geometric `L_loop` negligible | `L_J = 400-600 nH`, `L_loop ~ 20-50 pH`; the L_J sweep (5) is the check |
| 2 | quasi-static loop | loop ~25 um vs wavelength ~ mm-cm |
| 3 | PEC metal, no kinetic inductance (SQDMetal has no lambda_L) | `L_kin` folds into an effective `L_loop`; still << `2 L_J` |
| 4 | the integer in `Phi_loop = n Phi0` never enters | `c` is a fluctuation about the DC bias; the integer differentiates away |
| 5 | **NEW:** the array is a series chain of mutual inductances, `M_total = sum M_k` | true if the loops are galvanically in series and each threads the resonator flux independently — which is the SquidArray topology |

Imposed by hand (not from Palace): the step `phi = 2 pi Phi / Phi0`, the value
of `hbar` in `I_r^zpf`, and the junction nonlinearity. Palace supplies only
`f_r` and the modal current distribution (through `p`).

---

## 5. The `L_J` convergence check (analog of EPR_C_EXTRACTION.md 6.1)

The reference sweeps `L_probe` and fits `c_meas = c_inf * L_probe/(L_probe + L_loop)`.
The array analog, with the loop being a `2 L_J` circuit in series with its own
geometric `L_loop`:

```
c_meas(L_J)  =  c_inf * (2 L_J) / (2 L_J + L_loop)
```

Linearize: `1/c = 1/c_inf + (L_loop / c_inf) * 1/(2 L_J)`. Two `L_J` values give
`c_inf` (the geometric answer to report) and `L_loop` (a free byproduct).

**If `c(400 nH) ~ c(600 nH)` to within mesh noise**, then `L_loop << 2 L_J`,
the fit is unnecessary, and `c ~ c_inf` directly. That is what we see so far:
5-6% spread at every mesh density -> `c ~ c_inf`.

---

## 6. Mandatory validation (EPR_C_EXTRACTION.md 6)

- **6.1 `L_J` two-point** — done, spread 5-6% -> `c ~ c_inf`, no screening.
- **6.2 junction-field mesh convergence** — the B-field method (§2b) sidesteps
  the port-participation floor. `c_array` is stable to ~±5% over flux-cut
  heights 2-3 um; ~15% mesh sensitivity (box vs a1.5), finer = larger, not
  fully plateaued. See 7.4.
- **6.3 mode ID by participation, not frequency** — done. The resonator is the
  lowest-frequency mode with `sum_j |p_j| << 1` (array "phonon" modes have
  `sum_j |p_j| ~ 1`; a CPW or box mode has `<< 1e-3`). See
  `identify_resonator_mode` in `build_design.py`.
- `c < 0.1` sanity — per loop `c_k ~ 2e-3` (single-SQUID scale); `c_array ~
  0.05` is 23 healthy couplings summed coherently.
- **method self-check** — `M_k^loc` flat to ~11% across all loops, all `n`, all
  flux-cut heights.

---

## 7. Open questions

### 7.1 Which loops contribute? -- all of them for the loop-threading flux; the ENDS get extra lead coupling

`geometry_coupling.png`: the array runs down the exact centreline of a hairpin
fold of the resonator -- one leg 50 um above, the other 50 um below, both
full-length, joined by a U-turn at the far end. It also connects galvanically
to the launch pads at BOTH ends (`SquidArray` pins `a`, `b`).

**Loop-threading flux (B-field, `extract_M_bfield.py`):** `Phi_k` is a smooth
monotonic ramp, loop 0 (far end) `2.0e-13` -> loop 22 (resonator end)
`1.14e-13` Wb -- a factor ~1.8, not ~700x. `M_k = Phi_k/I_r ~ 1 pH` is FLAT
(0.98-1.72 pH). `c_k` runs `1.4e-4` -> `8.2e-5`, all same sign.
`c_array = sum_k c_k ~ 2.3e-3 rad` for n=23, `~ n * 9e-5` (linear). The user's
ParaView read -- "the eigenmode doesn't particularly have more flux on the
ends" -- is correct for this flux.

**The end loops -- current, not flux (`port-I.csv` + `port-EPR.csv`):**

| quantity | interior loops | end loops (0, 21, 22) |
|---|---|---|
| `M_k` from B-field flux | ~1.1 pH | 0.98-1.72 pH |
| `M_k` from `sqrt(2 L_J L_r p_EPR)` | ~1.8 pH (agree ~1.5x) | 11-46 pH (disagree ~30x) |
| junction current (`port-I.csv`) | ~1.5e-7 A | 1-4e-6 A |
| within-loop junction-current asymmetry | 0.3-0.5 | 0.03-0.13 |

The flux-derived and port-EPR-derived `M_k` now AGREE for the interior loops
(~1.5x, both ~1 pH). For the end loops port-EPR is ~30x higher. `port-I.csv`
shows why: the terminal-loop junctions carry ~10x more current, and their two
junctions carry EQUAL current (asymmetry -> 0) whereas interior loops have a
30-50% within-loop mismatch. The resonator-induced current circulates cleanly
through the chain-terminal loops (nowhere else to go) but partly BYPASSES the
interior junctions via the link bars / continuous array body. That extra
terminal current is inductive-energy participation, not loop-threading flux --
it does not change the geometric `M_k`.

**Earlier reading (retracted):** "only the two ends couple -- the array screens
its own interior." Wrong: the flux `M_k` is flat across all loops. port-EPR's
end weighting is a current-distribution effect at the chain terminals.

**Which `c`:** the flux-coupling `M_k` (-> `c ~ 2.3e-3`) is what modulates each
SQUID's phase in `-2 E_J cos(phi_dc/2 + c_k X) cos(...)`. The port-EPR total
(~4.7e-3) is closer to the array's inductive loading of the resonator.
Bracket: `c ~ 2-5e-3 rad`. See C_ARRAY_RESULTS.md §3-4.

### 7.2 The junction-corner field spike (RESOLVED -- does not matter)

ParaView close-up: `B_z` runs 2-3 orders of magnitude above the surrounding
loop-gap field at a point in many loops. It is a **sharp-corner (PEC)
singularity** at the 90-degree edges where the junction-port surfaces meet the
loop plates (`y ~ +/- 8-15 um`, on the loop-cell boundaries), a near-field
artefact that smooths away by `z ~ 2.5 um`.

The flux integral avoids it: the window (`|y| < 12.5 um`, inset from the
plates) excludes the corners. Dropping the brightest 1% of cells inside the
window changes `Phi_k` by < 2% at every height; max/median field ratio in the
window is only 1.2-2.4. Checked for loops 3, 11, 20 (`render_bfield_maps.py`).

### 7.3 Meissner screening vs the PEC model

Palace enforces `B.n = 0` on every metal surface and carries the surface
current `J_s` -- that IS the Meissner effect at the level of perfect screening.
The array body + the 3 um-away ground plane screen the resonator flux and the
result includes that (`M_k ~ 1 pH`, below the bare geometric value).

Not modelled: the London depth `lambda_L ~ 50-100 nm` / kinetic inductance
(SQDMetal has no `lambda_L`). It folds into an effective `L_loop` added to
`2 L_J`; with `L_J = 600 nH` and `L_kin ~ pH` it is a `< 1e-5` correction to
`c`. The `L_J` = 400 vs 600 check (§5, `c` flat to ~6%) bounds it. DC Meissner
/ trapped flux does not enter -- `c` is the small-signal fluctuation about the
DC bias and the integer in `Phi = n Phi0` differentiates away (assumption 4).

### 7.4 Absolute-value uncertainty

| source | affects | size |
|---|---|---|
| eigenmode field normalisation | -- | none (verified <1e-6 over 9 orders of s) |
| end-loop junction current (port-EPR sees it, B-field flux integral does not) | c (upward) | factor ~2 on c_array (2.3 -> 4.7e-3) |
| CPW current model (exact isolated-CPW vs twin-lead) | M only (cancels in c_k) | ~4x on M in henries |
| mesh density (box 2 um vs a1.5, 3.77 M elem) | c, M | 1.92e-3 -> 2.28e-3 (~19%) |
| flux-cut height (2-3 um) | c, M | ~3% |

`c_array ~ 2.3e-3 rad` (B-field flux coupling, all loops), bracket `2-5e-3`.
The shape (distributed, linear in n, flat `M_k`) is solid; the absolute
prefactor is good to a factor ~2.

### 7.5 Is `2 L_J` the right effective inductance?

The circulating (differential) loop mode sees `L_J + L_J = 2 L_J`; the common
mode sees `L_J / 2`. Flux coupling drives the circulating mode, so `2 L_J`.
(This enters the port-EPR path in §3; the B-field path in §2b does not use it.)

### 7.6 What does a single "c" mean for a 23-loop device?

The reference `c` is for one loop in one Hamiltonian term. Here each loop has
its own `c_k` in `-2 E_J cos(phi_dc/2 + c_k X) cos(phi_T^k)`. `c_array = sum
c_k` is the right scalar if the array is driven as one collective flux coupler;
for other operating points the per-loop `c_k` profile (§7.1) is the relevant
object.

---

## 8. What the code does

| step | file | output |
|---|---|---|
| build geometry, N loops, L_J per junction | `build_design.py` | qiskit-metal design |
| eigen solve, 2*N lumped inductive ports, pick resonator by `sum_j|p_j|` | `build_and_eigenmode.py` | `epr_result.json` with `p_per_junction` (2N) and `p_per_loop` (N) |
| **B-field: Phi_k, CPW-current inversion, M_k, c_k** | **`extract_M_bfield.py`** | **per-loop table + `.npz` raw arrays (the primary result)** |
| |B| / B_z colormaps over a loop + singularity-vs-height check | `render_bfield_maps.py` | PNGs |
| port-EPR path: `c_k`, `c_array`, noise filter, L_J fit, N-trend | `extract_c_array.py` | derived JSON (independent cross-check) |
| f + c convergence across AMR iterations | `check_convergence.py` | table + verdict |
| table + plots over all runs | `summarize_c.py` | `c_summary.json`, plots |
