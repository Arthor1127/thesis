# Deriving `c` from an EPR simulation, and why it depends on `L_probe`

This document derives, from first principles, (1) how the dimensionless flux-coupling constant
`c` follows from a Palace energy-participation-ratio (EPR) simulation, and (2) why the measured
value depends on the probe inductor `L_probe` used to make the measurement, culminating in the
exact functional form `c_meas(L_probe) = c_∞ · L_probe/(L_probe + L_loop)` that
`c_lprobe_fit.png` confirms to better than 0.02% over a full decade of `L_probe`.

---

## 1. The physical coupling: what `c` is

The qubit's SQUID loop is threaded by the readout resonator's current. The two junction phases
`φ_a`, `φ_b` combine in the circuit Hamiltonian as

```
H_JJ = -E_J cos(φ_a) - E_J cos(φ_b) = -2E_J cos((φ_a+φ_b)/2) · cos((φ_a-φ_b)/2)
```

Fluxoid quantization pins the *differential* phase to the flux threading the loop,
`φ_a - φ_b = 2π Φ_loop/Φ₀`. Write that flux as a DC bias plus the resonator's own zero-point
quantum fluctuation:

```
Φ_loop = Φ_dc + M·Î_r ,        Î_r = I_r^zpf (â_r + â_r†)
```

with `M` the mutual inductance between the resonator mode and the SQUID loop. Substituting,

```
(φ_a - φ_b)/2 = φ_dc/2 + c·(â_r + â_r†) ,      c ≡ π M I_r^zpf / Φ₀ = M I_r^zpf / (2φ₀)
```

`c` is therefore **purely a linear-circuit quantity** — a mutual inductance times a zero-point
current, divided by the flux quantum. It has nothing to do with the junction nonlinearity; it is
exactly the kind of number a classical Maxwell solver like Palace can compute, given the right
circuit boundary conditions. That is the entire justification for using an eigenmode field solve
instead of, say, a magnetostatic solve (which gives the wrong, DC current distribution) or a
lumped-capacitance model (which cannot represent the resonator's spatial current profile at all).

---

## 2. What Palace actually gives you: the participation ratio

You cannot measure `M` directly by leaving the SQUID's junction gaps open — with no closed current
path around the loop, no current flows and the "coupling" measured is identically zero. The
standard trick (Minev et al. 2021) is to **replace the junction with a known, finite lumped
inductor** `L_probe` at the gap (bridging the *other* junction, `jj_b`, with continuous metal so the
loop is a single-inductor circuit), and let Palace report how much of each eigenmode's inductive
energy ends up stored in that lumped port:

```
p_probe = (energy stored in the L_probe lumped port) / (total inductive energy of the mode)
```

This is a purely classical, well-defined, non-negative quantity that Palace computes natively via
its `LumpedPort` boundary condition (a surface across the junction gap relating tangential E-field
to current via `Z = R + iωL + 1/iωC`; with `R=C=0` this is an ideal inductor). It is written
directly to `port-EPR.csv`.

### 2.1 From participation to `c`

The port's stored energy is `(1/2) L_probe I_probe²`, where `I_probe` is the current flowing through
the loop (the current is the same everywhere in the loop, since it's a series circuit — see §4). By
definition of participation,

```
p_probe = (1/2) L_probe I_probe²  /  U_mode
```

Meanwhile the resonator mode's own zero-point current is fixed by its total energy,
`U_mode ↔ (1/2) ħω_m` quantum mechanically, or classically `U_mode = (1/2) L_r I_r²` for the
mode's lumped-equivalent inductance. Equating the *loop's* stored energy (built from `M` and `I_r`)
to the participation-derived energy, and using `E_probe/h ≡ φ₀²/(L_probe h)` as a convenient energy
scale for the probe inductor, algebra collapses to the extraction formula used throughout this
project:

```
c_m = (1/2) · sqrt( p_{m,probe} · f_m / (2 · E_probe/h) ) ,     E_probe/h = φ₀² / (L_probe · h) = 163.4 GHz / L_probe[nH]
```

(The leading factor of `1/2` traces back to the `(φ_a-φ_b)/2` in the trig identity in §1 — the
probe port's own zero-point phase swing is `2c`, not `c`.) This is the formula implemented in
`extract_epr_constant.py::c_from_participation`.

**Implicit assumption**: this formula is exact **only if the entire loop's inductance is the probe
inductor** — i.e. if the geometric metal trace forming the SQUID loop has negligible self-inductance
of its own. If that assumption fails, `p_probe` under-reports the loop's true participation, and
`c_m` computed this way underestimates the physical `c`. That is exactly what §4 quantifies.

---

## 3. Where the resonator's zero-point current comes from

For completeness (not needed to run the extraction, only to interpret it): if you want to recover
`M` itself as a cross-check,

```
M = 2 c φ₀ / I_r^zpf ,      I_r^zpf = sqrt( ħ ω_r / (2 L_r) )
```

with `L_r = 2 L_ℓ ℓ/π²` for a λ/4 CPW resonator of per-length inductance `L_ℓ` and length `ℓ`. None
of this enters the EPR extraction itself — Palace already accounts for the full spatial current
distribution automatically; this is only useful for sanity-checking the extracted `c` against a
hand analytic estimate.

---

## 4. Why `c_meas` depends on `L_probe`: the loop-inductance screening derivation

Model the SQUID loop as a series circuit: the geometric metal trace's own self-inductance,
`L_loop` (a property of the arm's length/width/gap — the quantity `compare_width_inductance.py`
estimates analytically), in series with the lumped probe inductor `L_probe` inserted at the junction
gap. The resonator's current `I_r` couples into this loop via mutual inductance `M`, driving an EMF
`ε = iωM I_r` around it. Since the loop's total size is tiny compared to the resonator's wavelength
at these frequencies (quasi-static regime — confirmed by the SQUID arm being at most ~10% of a
quarter-wavelength even at 600µm/arm), the loop behaves as a single lumped series `RLC`-like circuit
with impedance dominated by `iω(L_probe + L_loop)`:

```
I_loop = ε / (iω(L_probe + L_loop)) = M I_r / (L_probe + L_loop)
```

Crucially, **`I_loop` flows through both `L_probe` and `L_loop` equally** — same series circuit,
same current. The port's own stored energy is

```
E_probe = (1/2) L_probe I_loop² = (1/2) L_probe M² I_r² / (L_probe + L_loop)²
```

while the mode's total energy `U_mode ≈ (1/2) L_r I_r²` is (to leading, perturbative order)
unaffected by the small loop. So the participation ratio Palace reports is

```
p_probe = E_probe / U_mode  =  (M² / L_r) · L_probe / (L_probe + L_loop)²
```

**This is the key result: `p_probe` scales as `L_probe/(L_probe+L_loop)²`, not simply proportional
to `L_probe` alone.** As `L_probe → ∞` with `L_loop` fixed, the geometric loop inductance becomes
negligible next to the probe and `p_probe → M²/(L_r L_probe)`, recovering the naive
"no-screening" behavior; as `L_probe → 0`, essentially none of the loop's actual inductive energy
is captured by the tiny probe, and `p_probe → 0` even though the physical coupling `M` hasn't
changed at all.

### 4.1 Why `c_m` itself ends up *linear* in the screening factor

Now substitute this `p_probe` into the `c_m` formula from §2.1. Since
`1/(E_probe/h) = L_probe·h/φ₀² ∝ L_probe`, and `f_m` is essentially independent of `L_probe` (the
resonator frequency barely shifts for a weak, perturbative coupling),

```
c_m  ∝  sqrt( p_probe · L_probe )
     ∝  sqrt( [L_probe/(L_probe+L_loop)²] · L_probe )
     =  sqrt( L_probe² / (L_probe+L_loop)² )
     =  L_probe / (L_probe + L_loop)
```

The two `L_probe`-dependencies — `p_probe`'s intrinsic `L_probe/(L_probe+L_loop)²` screening, and the
extraction formula's own `1/(E_probe/h) ∝ L_probe` normalization (which exists precisely to convert
an energy fraction into a *phase* amplitude, i.e. to undo one power of the port's own energy scale)
— **combine to cancel one factor of `(L_probe+L_loop)`**, leaving exactly the observed law:

```
c_meas(L_probe) = c_∞ · L_probe / (L_probe + L_loop) ,     c_∞ ≡ M I_r^zpf/(2φ₀)  (the true, L_probe-independent answer)
```

This is not a curve-fit convenience — it is the *exact* consequence of modeling the loop as one
inductor (`L_probe`) in series with another (`L_loop`), both carrying the same current, one of which
you're measuring the energy in.

### 4.2 What this means physically

- **`L_loop ≪ L_probe`** (the regime you want when you *only* care about `c`): `c_meas ≈ c_∞`
  directly, no fit needed — this held for the original `l1=l2=300µm` geometry only marginally
  (`L_loop` 535–677 pH vs. `L_probe` 10–20 nH tested, giving 2.5–3.2% spread, passing the
  `EPR_C_EXTRACTION.md` §6.1 "negligible" check).
- **`L_loop` comparable to `L_probe`** (the `l1=l2=600µm` geometry, `L_loop ≈ 1.17 nH`): a *single*
  `c_meas` at any one `L_probe` is a biased underestimate of `c_∞` by the exact factor
  `L_probe/(L_probe+L_loop)`. The two/three-point fit isn't just a noise-averaging convenience here —
  it is *necessary* to recover the physical `c_∞`.
- **`L_loop` is a genuine, physically meaningful byproduct** of the fit, not a nuisance parameter:
  it is the SQUID loop's own self-inductance, extractable with no extra simulation cost. It roughly
  doubled (535–677 pH → 1170 pH) when the arm length doubled, consistent with the arm's
  length-dependent self-inductance derived in `lumped.ipynb`'s `line_inductance_isolated` /
  `line_inductance_cpw` formulas (both scale close to linearly with length in this regime).
- **Why `L_loop` in the "hundreds of pH to low nH" range is flagged** (`EPR_C_EXTRACTION.md` §8):
  SQDMetal's metals are PEC, so `L_loop` as fitted here is purely geometric — it omits the
  superconducting film's kinetic inductance. Once the fitted `L_loop` gets large enough to be a
  non-trivial fraction of `L_probe`'s design-relevant range (tens to hundreds of nH, per the typical
  transmon-junction-inductance scaling discussed separately), the *missing* kinetic contribution
  could itself matter for `β_L` and related screening physics, even though it does not change the
  validity of the `c_∞` extraction above (kinetic inductance would just be an additional term in
  series in the same derivation, folding into an effective `L_loop`).

---

## 5. Empirical confirmation

Fitting the `l1=l2=600µm` geometry's mode-1 data (`L_probe` = 10, 15, 20, 30, 50, 75, 100 nH)
to `c_meas = c_∞ L_probe/(L_probe+L_loop)` gives `c_∞ = 2.229×10⁻⁴`, `L_loop = 1172 pH`, with every
point matching the fitted curve to **better than 0.02%** (see `c_lprobe_fit.png`) — including points
spanning a full decade of `L_probe`, not just the two/three points nearest each other used for the
original fit. This is strong direct evidence that §4's derivation is the actual mechanism at play,
not a coincidental curve fit: the functional *form* itself is confirmed, not just its parameters.

---

## 6. Reference

- `EPR_C_EXTRACTION.md` (this directory) — the full extraction recipe, port setup, and validation
  checklist this document explains the physics behind.
- `C_EXTRACTION_RESULTS.md` (this directory) — the simulation log and numerical results this
  derivation is validated against.
- Minev et al., *Energy-participation quantization of Josephson circuits*, npj Quantum Information
  **7**, 131 (2021), arXiv:2010.00620.
