# ANSYS HFSS cross-check of the CPW impedance discrepancy — results

Independent third-solver check requested in `../HANDOFF.md`. Bottom line:

> **HFSS does not reproduce PALACE's ~40 Ω.** Its wave-port field impedance and
> its own S11-derived `Zin` for the bare CPW both land in the **50–55 Ω** band,
> right next to the two analytic formulas (≈50 Ω) and the PALACE
> capacitance-per-length result (47.6 Ω). PALACE's driven lumped-port
> `Zin = 40.5 Ω` (−19 %) is the outlier, exactly as the handoff hypothesised.
> The pattern repeats with the launchpad tapers (HFSS 44 Ω / PALACE 30 Ω) and
> the full design (HFSS 48 Ω / PALACE 21 Ω). Separately, an HFSS **eigenmode**
> solve puts the resonator at **f_r = 8.55 GHz**, ~3 % above PALACE's own
> eigenmode value — i.e. the resonance is where PALACE says it is.

## Method (all tests)

- ANSYS Electronics Desktop **2024 R2**, HFSS **Driven-Modal**, driven via
  `ansys-aedt-core` 1.5 (COM transport — gRPC does not connect on this machine).
- Geometry rebuilt from the handoff GDS metal layers (datatype 0 = ground,
  10 = pads, 11 = trace); the cheese-hole datatypes are ignored. The ground
  plane's "keyhole" polygons are healed to proper polygons-with-holes
  (`shapely.buffer(0)`) before being drawn.
- **Metal: zero-thickness Perfect-E** sheets (matches SQDMetal's quasi-3D
  treatment, as requested).
- Silicon slab `εr = 11.45`, loss-free; **grounded Perfect-E backplane 500 µm
  below the trace** (matches PALACE `z_prop_bottom = 0` and the
  grounded-backplane analytic formula).
- **Wave ports** on the two ±x domain end-faces, spanning the full
  substrate+vacuum cross-section, with an **integration line across the CPW
  gap** so HFSS reports the field-based characteristic impedance `Zo` (= `Zpi`).
- S-parameters exported **non-renormalised** (reference = each port's `Zpi(f)`),
  then renormalised in post to a **flat 50 Ω** so that
  `Zin = 50·(1+S11)/(1-S11)` is computed exactly the way PALACE computes it
  (`Z0_ref = 50`).
- Adaptive mesh to `ΔS ≤ 0.02` (2 consecutive); interpolating sweep 6–10 GHz.

Reproduce: `python run_ansys_tests.py A B C` in the `qiskit_qc` env.

---

## Test A — bare CPW line (highest priority)

`geometry/line_only.gds`: a 2 mm straight CPW (w = 20 µm, gap = 11.5 µm), open
ends, no launchers. Rebuilt as clean rectangles; ground-plane half-width swept
150 → 1000 µm to check domain sensitivity.

| ground half-width | HFSS field **Zpi** (port 1 / port 2) | **Zin from S11** (renorm. 50 Ω) | \|S11\| @ 50 Ω | ε_eff | converged |
|---|---|---|---|---|---|
| 150 µm  | 54.7 / 54.8 Ω | 51.9 + 2.0j Ω | 0.027 (−31 dB) | 6.224 | pass 11, ΔS 0.013 |
| 200 µm  | 54.2 / 56.4 Ω | **49.6 + 1.4j Ω** | 0.015 (−37 dB) | 6.224 | pass 10, ΔS 0.015 |
| 400 µm  | 55.1 / 56.3 Ω | 52.0 + 3.5j Ω | 0.040 (−28 dB) | 6.200 | pass 11, ΔS 0.012 |
| 1000 µm | 57.9 / 56.6 Ω | 56.5 + 6.1j Ω | 0.084 (−21 dB) | 6.270 | pass 10, ΔS 0.015 |

Observations:

- **ε_eff ≈ 6.22 = (εr + 1)/2 exactly** in every run — the textbook
  thick-substrate CPW value. HFSS is modelling the cross-section correctly.
- **Zpi is flat vs frequency** (std < 0.05 Ω across 6–10 GHz) and **stable at
  ~54–56 Ω** across a 6× range of ground width — the ~10 % offset from the
  50.3 Ω analytic value is *not* a finite-ground or domain-truncation artifact.
  (It is within the known spread of CPW `Zpi` vs the specific analytic
  conformal-mapping definition; the wider-domain 1000 µm point drifts up a
  little and is the least trustworthy — larger air box, same pass count.)
- **`Zin` from HFSS's own driven S11, computed identically to PALACE**, is
  **≈ 50–52 Ω** for the compact domains, and the raw \|S11\| renormalised to
  50 Ω is **−28 to −37 dB** — i.e. HFSS's driven simulation sees a CPW line
  that is *well matched to 50 Ω*.
- PALACE's driven sim, same geometry, same `Zin` formula, gave
  **40.45 − 1.66j Ω** with \|S11\| far larger.

### How the methods line up (bare cross-section)

| method | Z0 / \|Zin\| | deviation from 50 Ω |
|---|---|---|
| DesignHelper conformal-map | 49.96 Ω | −0.1 % |
| SQDMetal `CpwParams` (grounded backplane) | 50.29 Ω | +0.6 % |
| PALACE capacitance-per-length | 47.6 Ω | −4.9 % |
| **HFSS `Zin` from driven S11 (50 Ω ref)** | **≈ 50–52 Ω** | **±3 %** |
| **HFSS wave-port field `Zpi`** | **≈ 54–56 Ω** | **+8–11 %** |
| **PALACE driven lumped-port `Zin`** | **40.5 Ω** | **−19 %** |

Every independent estimate except PALACE's driven lumped-port `Zin` sits
within ~10 % of 50 Ω. **PALACE's driven-port `Zin` is confirmed as the
outlier.**

Files: `testA_line_only.s2p` (raw, ref = Zpi), `testA_line_only_50ohm.s2p`
(renormalised), `testA_line_only_sparams.csv`, `testA_line_only_Zo.csv`,
`testA_line_only_convergence.csv`; likewise `testA_yh150/400/1000*`.

---

## Test B — CPW + launchpad tapers (no coupler)

`geometry/with_launchpads.gds`: same CPW trace with a real
launchpad/wire-bond-pad taper at each end; wave ports on the **outer bond-pad
edge** (x = ±1750 µm).

| quantity | HFSS | PALACE |
|---|---|---|
| field `Zpi` (port 1 / 2, at the pad) | 45.0 / 45.2 Ω | — |
| **`Zin` from S11 (50 Ω ref), 6–10 GHz mean** | **44.1 − 0.2j Ω** | **≈ 29.7 Ω** |
| \|S11\| @ 50 Ω | 0.062 (−24 dB) | (much larger) |
| ε_eff (pad port mode) | 5.94 | — |
| converged | pass 13, ΔS 0.012, 27 k tets | — |

- Both solvers agree the taper **pulls the apparent `Zin` down** from the bare
  line.
- HFSS: 50 → 44 Ω, i.e. the taper contributes about **−6 Ω (−12 %)**.
- PALACE: 50 → 29.7 Ω, i.e. the taper contributes about **−20 Ω (−40 %)** in
  their sim.
- So HFSS sees the taper as a real but modest contributor, not the dominant
  source of mismatch. The gap between the two solvers *widens* with the taper
  added, consistent with the PALACE driven lumped-port method reading
  systematically low and the error compounding.

Files: `testB_launchpads*`.

---

## Test C — full hanger-resonator design

`geometry/full_design_cs14_cl300.gds`, 2-port driven S21. The feedline
(+ bond pads) is the driven signal conductor; the hanging resonator meander is
placed in the ground net (it is short-circuited to ground through its
`ShortToGround` and only capacitively coupled to the feedline). Wave ports on
the outer bond-pad edges (x = ±1750 µm). The resonator hangs to y ≈ −1317 µm,
so the domain has to be **y_half ≈ 1570 µm** to contain it (an earlier run
that clipped at ±600 µm chopped the resonator and it did not resonate).

### Baseline (away from resonance) — robust

| quantity | HFSS | PALACE |
|---|---|---|
| field `Zpi` (port 1 / 2) | 46.0 / 45.9 Ω | — |
| **`Zin` from S11 (50 Ω ref), band mean** | **47.0 − 0.3j Ω** | **≈ 20.9 + 5.5j Ω** |
| \|S11\| @ 50 Ω | 0.035 (−29 dB) | (much larger) |
| ε_eff | 6.24 | — |
| converged | pass 9, ΔS 0.012, 69 k tets | — |

The off-resonance `Zin` of the full design is **47 Ω in HFSS vs 20.9 Ω in
PALACE** — the same pattern as Tests A and B, now with the largest absolute
gap. Three geometries, three confirmations that HFSS's driven simulation does
not reproduce PALACE's low `Zin`.

### The resonance — from an eigenmode solve

**`f_r` (HFSS eigenmode) = 8.546 GHz** vs PALACE's 8.101 GHz (fitted notch) /
8.25–8.28 GHz (their own eigenmode). HFSS is ~3 % above PALACE's eigenmode,
~5 % above their fit — a normal solver/mesh-level difference for a distributed
resonator, and it confirms the resonator sits where PALACE says it does.

| eigenmode | HFSS | identification |
|---|---|---|
| mode 1 | **8.546 GHz**, Q ≈ 92 000 | the hanging resonator |
| modes 2–5 | 16.2 / 16.9 / 18.7 / 20.5 GHz | feedline + box modes |

`Q ≈ 92 000` for mode 1 equals `1/tan δ` (with `tan δ = 1.2e-5`), i.e. the
resonator is dielectric-loss-limited exactly as set up — a good consistency
check.

**`Q_c` (the coupling Q) was not extracted.** The driven S-matrix is flat: the
**raw** `|S11|` stays `0.054 ± 0.002` across 6.5–8.7 GHz with no resonance
feature, and every driven method agrees (interpolating sweep, Fast/ALPS ROM,
fine discrete to 80 kHz, and a solve adapted *at* `f_r`). The driven adaptive
mesh converges on the smooth through-line transmission and never resolves the
coupling capacitance, which is small: the `CoupledLineTee` puts a **~15 µm
strip of ground metal between the feedline and the resonator's coupled arm**
(`testC_coupler_geometry.png`), so the two couple only weakly across
gap + ground + gap.

Renormalising the raw S-parameters to 50 Ω through the HFSS "Port Zo" report
(which spikes to 2 Ω / 100 Ω at a single frequency near the mode) *injects*
fake `S21` notches — those are artefacts, not physics; after de-glitching, only
a ~0.7 dB residual blip at 8.43 GHz remains, a faint corroboration of a
resonance near the eigenmode value.

`Q_c` does not affect the impedance conclusion, and the `f_r` cross-check above
is the meaningful resonance result. (The handoff itself notes this design shows
"a mild double-dip" in PALACE's ported driven-S21 — consistent with a
resonator that couples weakly and interacts with feedline/box modes.)

Files: `testC_eigenmodes.json`, `testC_eigen.py`, `testC_full_design.s2p`
(driven, adapted at `f_r`), `testC_coupler_geometry.png`.

---

## Deviations from the handoff (all deliberate)

1. **Compact domain box** — a few hundred µm of silicon / vacuum / ground
   around the structure, not the full 4.8 × 2 mm chip. Justified by the
   handoff's own finding that lateral ground width (0.1×–3×) and vacuum height
   (1×–12×) have zero effect; a compact uniform cross-section is the standard
   choice for wave-port field accuracy. Test A's ground-width sweep re-confirms
   the insensitivity.
2. **No wire-bonds / air-bridges** modelled (planar only). The wave port
   excites the CPW mode symmetrically, so the slot-line mode the bonds would
   suppress is not excited.
3. **Zero-thickness Perfect-E metal** (as requested).
4. **Wave ports on the domain end-faces**, so the modelled ground ends at the
   port plane rather than extending past it. Test A shows this does not matter.
5. HFSS reports `Zpi` (integration-line power–current definition); `Zpv` / `Zvi`
   were not separately extracted (they require re-solving with a different port
   setting). For a quasi-TEM line the three differ by only a few percent, far
   less than the PALACE discrepancy.
