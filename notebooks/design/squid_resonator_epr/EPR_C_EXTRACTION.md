# Extracting the flux-coupling constant `c` from a Palace EPR eigenmode simulation

**Audience:** an agent with access to the design model and the SQDMetal/Palace stack, but no
prior context on why this quantity is wanted.

**Deliverable:** a single dimensionless number `c` (radians), per resonator mode, plus an
uncertainty and a set of passed/failed validation checks.

---

## 1. What `c` is and why anyone cares

The device is a transmon qubit whose junction connects to ground, coupled to a CPW resonator
**not capacitively but through a large symmetric SQUID loop**. The resonator's current threads
magnetic flux through that loop. That is the entire coupling mechanism.

The two SQUID junctions contribute a term to the circuit Hamiltonian that factorizes exactly
by trigonometric identity:

```
-E_J cos(φ_a) - E_J cos(φ_b)  =  -2 E_J cos((φ_a+φ_b)/2) cos((φ_a-φ_b)/2)
```

Fluxoid quantization ties the differential phase to the loop flux, `φ_a - φ_b = 2π Φ_loop/Φ₀`.
Splitting the loop flux into a DC bias plus the resonator's quantum fluctuation gives

```
(φ_a - φ_b)/2  =  φ_dc/2  +  c·X̂ ,        X̂ = â_r + â_r†
```

so the Hamiltonian term becomes `-2 E_J cos(φ_dc/2 + c·X̂) cos(φ̂_T)`, and

```
c = π M I_r^zpf / Φ₀ = M I_r^zpf / (2 φ₀)
```

with `M` the mutual inductance between the resonator mode and the SQUID loop, `I_r^zpf` the
resonator's zero-point current, and `φ₀ = Φ₀/2π = 3.29106 × 10⁻¹⁶ Wb`.

`c` is the **only** quantity in that term that requires a field solver. `E_J`, `φ_dc`, and the
transmon's `φ_T^zpf` are all set by design or by the experimentalist. Everything downstream —
notably the three-wave mixing rate `g₃ = E_J sin(φ_dc/2) · c · (φ_T^zpf)²` — is algebra once
`c` is known.

**Expected magnitude:** `c` is small, of order 10⁻³ to 10⁻² rad. A value above ~0.1 means
something is wrong (see §7).

---

## 2. Why this is done with EPR and not a magnetostatic solve

Reasonable-sounding alternatives that do **not** work:

| Approach | Why it fails |
|---|---|
| Magnetostatic solve for `M` | Gives a DC uniform current distribution. A λ/4 resonator has a cosine current profile; the DC `M` is not the modal `M`. |
| Capacitance matrix / LOM | Collapses the resonator to one lumped node. Cannot represent position along the line, cannot represent harmonics. |
| Remove the SQUID junctions, leave gaps open | No circulating current is possible. Participation is identically zero. Measures nothing. |
| Bridge both junction gaps with metal | The loop becomes a perfect superconducting screen and expels the resonator's flux. Measures the *screened* response, not the bare coupling. |

The working method replaces the junctions with a **finite lumped inductor** and reads the
energy-participation ratio, which Palace computes natively.

---

## 3. Physics assumptions — read this before trusting any number

Palace is a classical Maxwell solver. It has no ħ, no Φ₀, no order parameter. It does **not**
know the metal is a superconductor. The division of labour is:

- **From Palace (classical):** mode frequencies `ω_m`, and the current distribution that
  determines `M`. Flux conservation in a zero-resistance loop is pure Faraday and is handled
  correctly by a PEC loop with a lumped inductor.
- **Imposed by hand (quantum):** the step `φ = 2π Φ/Φ₀`, the value of `ħ` in `I_r^zpf`, and the
  junction nonlinearity.

The integer in `Φ_loop = nΦ₀` never enters, because `c` describes a *fluctuation* about the DC
bias, and the integer differentiates away. This is why a classical solver suffices.

**Known omission:** metals are PEC, so kinetic inductance from the film's London penetration
depth is not modelled. SQDMetal has no λ_L support at all (verified: zero occurrences in the
source tree; the `Material` class carries only permittivity, permeability, loss tangent, and
raises on unknown kwargs). Consequence: the extracted loop inductance is geometric only, hence
underestimated, hence `c` slightly overestimated. With `L_probe = 10 nH` and `L_loop ~ 50 pH`
the error on `c` is well under 1% — negligible. State it as an assumption; do not attempt to
patch it. It matters only for `β_L` (§8).

---

## 4. Simulation setup

### 4.1 Port assignment

| Port | Physical element | Value |
|---|---|---|
| `J_T` | transmon junction to ground | `L_T` (design value) |
| `J_probe` | **single** lumped inductor replacing the SQUID | `L_probe`, default 10 nH |

The SQUID has two junction gaps in the geometry (`jj_a`, `jj_b`). Configure them as:

- `jj_a` → lumped port, `R = 0`, `C = 0`, `L = L_probe`
- `jj_b` → **bridged with metal** (short)

Both arms must not be left open. One arm shorted plus one arm carrying a finite inductor makes
the loop a single-inductor circuit, which is the configuration where the port's branch flux
equals the loop flux.

**Do not attach a 50 Ω port to the feedline for these runs.** Resistive ports go to `port-Q.csv`,
not `port-EPR.csv`, and loading the resonator only degrades the eigenmode.

### 4.2 How a lumped port works mechanically

A lumped port in Palace is a **boundary condition on a surface**, not a circuit element dropped
into a netlist. The junction gap is a vacuum break in the PEC metal. A rectangle spanning that
gap is tagged as its own physical surface in the mesh; Palace relates the tangential E across
that surface to the current through it via `Z = R + iωL + 1/(iωC)`. With `R = C = 0` the gap
behaves as an ideal inductor bridging the two metal edges.

`Direction` must point **across** the gap, along the current path. A wrong direction yields a
port that couples to nothing.

### 4.3 Solver settings

- Solve **at least 6 modes** — the resonator fundamental plus harmonics, and the transmon mode.
- Choose `L_T` so the transmon mode sits **1–2 GHz away** from the resonator. Hybridization is
  not needed for this measurement and only muddies mode identification.
- Adaptive mesh refinement off for production runs unless separately validated.

---

## 5. The extraction

Palace writes `port-EPR.csv` containing `p_{m,j}` — the fraction of mode `m`'s inductive energy
stored in inductive port `j` — for every mode and every inductive lumped port.

For each resonator mode `m`:

```
c_m = (1/2) · sqrt( p_{m,probe} · f_m / (2 · E_probe/h) )
```

with everything in frequency units and

```
E_probe/h = φ₀² / (L_probe · h) = 163.4 GHz / L_probe[nH]
```

The factor of ½ traces back to the `(φ_a - φ_b)/2` in the trig identity — the probe port's
zero-point phase is `2c`, not `c`.

**Worked check.** At `f = 5.039 GHz`, `L_probe = 10 nH` (so `E_probe/h = 16.34 GHz`):

```
c = 0.196 · sqrt(p)
```

`p = 3×10⁻³ → c ≈ 0.011 rad`. `p = 10⁻⁴ → c ≈ 0.002 rad`. Use this to sanity-check any
implementation before trusting it.

**Convention note:** this `c` multiplies `(â + â†)`. If a downstream consumer defines
`φ_resonator` with its own zero-point factor already included, divide accordingly. The two
differ by a factor of order 0.05 — get this wrong and everything downstream is silently off.

### 5.1 Optional: recover `M`

Not needed for the Hamiltonian, useful only as a cross-check against analytics:

```
M = 2 c φ₀ / I_r^zpf ,     I_r^zpf = sqrt(ħ ω_r / (2 L_r))
```

`L_r` is the mode's lumped-equivalent inductance; for a λ/4 CPW, `L_r = 2 L_ℓ ℓ / π²`.

---

## 6. Required validation

### 6.1 Two-point `L_probe` fit — mandatory

Run the identical geometry at `L_probe = 10 nH` and `L_probe = 20 nH`. The measured value obeys

```
c_meas(L_probe) = c_∞ · L_probe / (L_probe + L_loop)
```

Two runs give two equations in two unknowns. Solve for:

- `c_∞` — the geometric answer, the one to report
- `L_loop` — the SQUID loop's own inductance, a free byproduct

If the two runs agree to within mesh noise, `L_loop ≪ L_probe` and `c ≈ c_meas` directly.
If they disagree substantially, the loop is screening and the fit is doing real work.

This check is what converts "I hope `c` is geometric" into "I measured that it is." Do not skip it.

### 6.2 Junction-port mesh convergence — mandatory

`c ∝ √p`, and `p` is the current through a ~3 µm gap in a multi-millimetre model. **The
eigensolver's convergence criteria say nothing about junction field accuracy.** A run can report
backward error 10⁻¹⁸ on a mesh where `p` is wrong by a factor of two.

Sweep mesh density local to the junction ports and plot `p` vs. element count. Report `c` only
once `p` has plateaued. Because of the square root, an order of magnitude in `p` is only a
factor of 3 in `c` — forgiving — but `p` still needs to be good to better than ~2× for `c` to be
good to 40%.

### 6.3 Mode identification by participation, not frequency

The transmon mode has `p_T ≈ 1`. Resonator modes have `p_T ~ 10⁻²`. **Sort on participation.**
Sorting by eigenfrequency fails silently the moment a qubit mode crosses a resonator harmonic,
and produces a plausible-looking wrong number.

Record which mode index maps to which physical mode explicitly in the output.

---

## 7. Failure modes and what they look like

| Symptom | Cause |
|---|---|
| `p_probe` exactly 0 | Both junction gaps left open; no current path around the loop |
| `p_probe` implausibly large, `c > 0.1` | `jj_b` not actually shorted, or `L_probe` far too small |
| `c` drifts strongly with `L_probe` | `L_loop` comparable to `L_probe`; use the §6.1 fit, consider larger `L_probe` |
| `c` drifts with mesh density | Not converged; §6.2 |
| No `port-EPR.csv` written | No **inductive** lumped ports defined, or only resistive ports present |
| Port present but `p` ~ 0 | `Direction` not across the gap |

---

## 8. Downstream consumers (context, not tasks)

- **Three-wave mixing:** `g₃ = E_J sin(φ_dc/2) · c · (φ_T^zpf)²`, with
  `φ_T^zpf = (2E_C/E_J,T)^(1/4) ≈ 0.3`.
- **Optimal bias:** linear coupling ∝ `sin(π Φ_dc/Φ₀)`, maximal at `Φ_dc = Φ₀/2`, where the
  static coupling simultaneously vanishes — pure parametric operation. In practice the bias sits
  near but not at `Φ₀/2`, limited by fabricated junction asymmetry (typically a few percent).
- **Screening parameter:** `β_L = 2π L_loop I_c / Φ₀`, using `L_loop` from the §6.1 fit. If
  `β_L` approaches 1 the flux-tuning curve becomes hysteretic and the small-`c` linearization
  is unsafe. **This is the one place where the missing kinetic inductance (§3) matters** — if
  the fit returns `L_loop` in the hundreds of pH rather than tens, escalate rather than
  proceeding.

---

## 9. Output contract

Keep the FEM run and the algebra in separate programs. The solve is expensive (~4000 s wall
clock at ~4×10⁶ unknowns, ~32 GB); the algebra runs in a second. Dump raw quantities only, and
re-derive `c` offline without re-solving.

```json
{
  "l_probe_nh": 10.0,
  "l_transmon_nh": 12.0,
  "modes": [
    {"m": 1, "f_ghz": 5.0389, "Q": 4.03e5,
     "p": {"probe": 3.1e-3, "transmon": 1.2e-2},
     "assignment": "resonator_fundamental"}
  ],
  "mesh": {"elements": 3912808, "min_size_um": 2.0, "junction_min_size_um": 0.5}
}
```

The derived file should carry `c_m` per mode, `c_∞`, `L_loop`, the residual of the two-point
fit, and an explicit pass/fail on each check in §6.

---

## 10. Reference

- Minev et al., *Energy-participation quantization of Josephson circuits*, npj Quantum
  Information **7**, 131 (2021), arXiv:2010.00620 — the method.
- SQDMetal paper arXiv:2511.01220, Appendix C — a worked transmon-resonator EPR example in this
  exact software stack, with numbers to reproduce as a validation case.
