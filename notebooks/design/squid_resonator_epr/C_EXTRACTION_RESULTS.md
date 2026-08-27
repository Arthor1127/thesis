# Flux-coupling constant `c`: extraction results and simulation log

**Date:** 2026-08-24/25
**Design:** `build_design.py` — circular transmon (`CircTransmonSQUID`) with an added SQUID loop,
coupled to a meandered CPW readout resonator (`readout_1`), on the local workstation and the CAB
cluster (`arturo.ruiz@10.73.25.223`, SGE).

---

## 1. What was measured, and why

The qubit couples to the resonator not capacitively but through a SQUID loop threaded by the
resonator's current. The two SQUID junctions' Hamiltonian term factorizes by trigonometric
identity, and fluxoid quantization ties the differential junction phase to the loop flux. Splitting
that flux into a DC bias plus the resonator's zero-point fluctuation defines a single dimensionless
constant `c`:

```
-E_J cos(φ_a) - E_J cos(φ_b) = -2 E_J cos((φ_a+φ_b)/2) cos((φ_a-φ_b)/2)
(φ_a - φ_b)/2 = φ_dc/2 + c·X̂,     X̂ = â_r + â_r†
c = π M I_r^zpf / Φ₀ = M I_r^zpf / (2 φ₀)
```

`c` is the only quantity in the coupling term that needs a field solver — `E_J`, `φ_dc`, and the
transmon's zero-point phase are set by design. It is extracted from Palace's energy-participation
ratio `p_probe` (fraction of a mode's inductive energy stored in the SQUID's lumped inductive port),
not from a magnetostatic mutual inductance (which would give the DC current distribution, not the
resonator mode's cosine-profile current) or a lumped capacitance model (which collapses the
resonator to one node and cannot represent its harmonics):

```
c_m = (1/2) · sqrt( p_{m,probe} · f_m / (2 · E_probe/h) ),     E_probe/h = φ₀²/(L_probe·h) = 163.4 GHz / L_probe[nH]
```

Full derivation, the port setup (one SQUID arm bridged with metal, the other a lumped inductor
`L_probe`), and the mandatory validation checks are in `EPR_C_EXTRACTION.md` §1–6; this document
only summarizes what was actually run against that recipe. Two checks from that spec were carried
out here: the §6.1 two-point (here three-point) `L_probe` fit, and mode identification by
participation rather than frequency (§6.3). The §6.2 junction-mesh-density convergence check was
partially addressed (see §5 below) but not completed as a formal sweep.

---

## 2. Simulation history

### 2.1 First working baseline (local, coarse)

A first coarse smoke test (`run_coarse_test.py`, uniform mesh: `readout_1` and the qubit component
each `min_size=25µm`, no per-region splitting, `solver_order=1`, `starting_freq=7.0GHz`) ran
successfully in about a minute, 434,131 elements, and gave the first working `c` values:

| Mode | f (GHz) | c |
|---|---|---|
| resonator_mode_1 | 7.748 | 2.00×10⁻⁴ |
| resonator_mode_2 | 14.563 | 1.85×10⁻³ |

This became the reference "known-good" mesh/solver configuration for everything that followed.

### 2.2 Finer local mesh attempt → OOM

Editing the smoke test to `readout min=4µm`, `qubit min=2.5µm` caused mesh generation to run for
20+ minutes and climb past 38 GB on a 14 GB-RAM workstation. Moved all further work to the cluster.

### 2.3 Per-region mesh-size sweep (cluster)

To find where GMSH's background-field computation becomes impractical, each of four regions
(`readout`, `pad`, `squid_arm`, `junction`) was swept independently (min_size 25→1µm, others held
at the 25µm baseline), each step capped at a 2-minute timeout:

| min_size | all 4 regions |
|---|---|
| 25–12µm | pass, 21–27s |
| 10–6µm | pass, 22–96s (element count climbing from ~700k to ~6M) |
| **5µm and below** | **timeout, every region, no exceptions** |

This is a shared GMSH background-field-computation floor around 6µm — not a property of any one
region's size or shape. Combining multiple regions below 6µm simultaneously (tested at "6µm
everywhere") produced a 5.96M-element mesh with **38.6M unknowns** (`solver_order=2`) and a
**119 GB** measured peak — confirming that stacking multiple sub-6µm regions compounds cost
sharply, which is what had caused the original local OOM.

### 2.4 Adaptive mesh refinement (AMR)

Rather than manually guessing which regions need finer resolution, Palace's own AMR
(`enable_mesh_refinement`, Dörfler marking, `SaveAdaptMesh=True`) was used to let the solver decide,
starting from the proven 25µm uniform mesh.

**First AMR run** (`solver_order=1`, `Target=7.0GHz`, `MaxIts=3`, conformal refinement —
non-conformal refinement segfaulted on this mesh and was not pursued further):

| Iteration | Unknowns (ND) | Indicator norm | Wall time | Peak memory |
|---|---|---|---|---|
| 1 | 499,472 | 3.334×10⁻¹ | 53.7s | 11.8 GB |
| 2 | 1,588,408 | 2.368×10⁻¹ | 207.5s | 38.2 GB |
| 3 | 6,154,427 | 1.712×10⁻¹ | 866.6s | 160.4 GB |
| 4 (attempted) | 23,021,876 | — | — | **>513 GB, node swap-thrashed** |

Iteration 4 was killed manually once memory exceeded the 503 GB physical node and it began
swapping; the node recovered fully and iterations 1–3's data were untouched.

**Where AMR actually refined.** Parsing the per-element `Indicator` field (Palace's own error
indicator, written to a dedicated ParaView cycle — see `analyze_indicator_field.py`) against
hand-defined region boxes showed the qubit pad and SQUID arms together accounted for only
~1–4% of total indicator mass, while the highest-indicator points clustered tightly along the
**readout resonator's trace edges**. Field plots of the two converged modes confirmed this
directly: mode 1 (7.97 GHz)'s |E| lives entirely on the resonator meander with essentially nothing
on the qubit; mode 2 (14.5 GHz)'s field is dominated by the **chip boundary walls** (a rectangular
frame pattern) — i.e. **mode 2 is a spurious box/cavity mode** of the finite simulation domain, not
a real second resonator harmonic, and its `c` value is not physically meaningful.

### 2.5 Mesh re-tuned from the AMR result, second AMR run

Based on §2.4, the manual mesh was re-tuned: `readout_1`'s fine-mesh taper widened
(`taper_dist_min` 10→40µm, giving the fine zone more reach along the resonator) and the qubit's
narrowed (20→8µm, since AMR showed it needed little). A quick direct check of the resulting mesh
file confirmed the SQUID loop was still resolved at ~3.1µm median element size from the qubit-wide
pass alone, with no dedicated region needed there.

Re-running AMR (`MaxIts=2`, same base recipe) on this re-tuned mesh:

| Iteration | Unknowns | Indicator norm | Peak memory |
|---|---|---|---|
| 1 | 2,863,284 | 2.699×10⁻¹ | 79.6 GB |
| 2 (converged, final) | 4,579,392 | 1.972×10⁻¹ | 129.2 GB |
| 3 (attempted) | — | — | crashed: `mfem::ErrorException` in `STable3D` (mesh-topology bug during refinement rebuild, not a memory issue) |

Iteration 2's converged result (safely within the 150 GB budget requested) gave:

| Mode | f (GHz) | Q | c |
|---|---|---|---|
| resonator_mode_1 | 7.968 | 2.65×10⁷ | 8.67×10⁻⁵ |
| resonator_mode_2 (box mode) | 14.522 | 4.79×10⁵ | 1.676×10⁻³ |

Mode 1's `c` dropped by more than 2× from the coarse-mesh value (2.00×10⁻⁴ → 8.67×10⁻⁵),
confirming the coarse mesh had been significantly under-resolving the readout resonator; mode 2's
`c` barely moved (1.85×10⁻³ → 1.68×10⁻³), consistent with it being governed by chip-boundary
geometry rather than local mesh density.

### 2.6 Mesh-blowup debugging (build_and_eigenmode.py)

Applying the same taper re-tuning to the production script (`build_and_eigenmode.py`) for a
manual (non-AMR) `L_probe` sweep hit three compounding mesh-size bugs before converging on a safe
configuration, isolated one at a time via controlled A/B mesh-only checks:

1. A dedicated `fine_mesh_in_rectangle` added over the SQUID-arm bounding box fine-meshed the
   *entire box* (mostly empty substrate/ground plane) regardless of where actual metal was —
   unlike `fine_mesh_components`, which only refines near real geometry edges. Pushed unknowns
   from 2.86M to 13.2M (peaked at 156 GB on a 188 GB node before being killed). **Removed** — the
   qubit-wide pass alone already resolves the loop adequately (§2.5).
2. `--qubit-bg-max-size-um` defaulted to 20µm (not 150µm), forcing a large fraction of the *whole*
   domain to stay at 20µm instead of relaxing to background coarseness. **Default raised to 150µm.**
3. The existing per-junction-port `fine_mesh_in_rectangle` (tiny boxes, `min=1.5µm` default) was,
   on its own, *also* below the §2.3 6µm floor and, combined with the other fine regions, pushed a
   otherwise-52.9s/2.47M-element mesh past a 3-minute timeout — confirmed by isolation testing
   (identical config with/without this block: 52.9s vs. timeout). **Removed entirely**; junction
   resolution is already covered by the qubit-wide pass.

The corrected recipe (readout `min=6µm, taper=40µm`; qubit `min=2.0µm, max=150µm, taper=5µm`; no
junction- or squid-specific rectangle) meshes in ~30–90s at ~1.0–2.5M elements depending on exact
settings, safely reproducible.

A separate, unrelated bug was also found and fixed: `build_and_eigenmode.py`'s `_json_safe()`
helper converted numpy complex scalars via `.tolist()`, which yields a plain Python `complex` —
still not JSON-serializable — corrupting `epr_result.json` mid-write on a run whose solve had
already succeeded. Fixed by recursively re-checking the `.tolist()` result.

### 2.7 `L_probe` sweep — the mandatory §6.1 check

With the corrected mesh (1,280,457 unknowns, ~31.5 GB peak per run, well inside budget), the design
was solved at three values of `L_probe` (10, 15, 20 nH), fitting `c_meas(L_probe) = c_∞·L_probe/(L_probe+L_loop)`:

| L_probe | c_meas, mode 1 (resonator) | c_meas, mode 2 (box mode) |
|---|---|---|
| 10 nH | 8.223×10⁻⁵ | 1.604×10⁻³ |
| 15 nH | 8.364×10⁻⁵ | 1.638×10⁻³ |
| 20 nH | 8.437×10⁻⁵ | 1.656×10⁻³ |

**Mode 1 fit:** `c_∞ = 8.66×10⁻⁵`, relative spread across the three points **2.5%** (passes the
"loop inductance negligible" check), `L_loop = 535 pH`. This matches the independent AMR-2 estimate
(8.67×10⁻⁵) almost exactly, despite coming from a different mesh and pipeline — strong convergence
evidence.

**Mode 2 fit:** `c_∞ = 1.712×10⁻³`, spread 3.2% — numerically stable too, but not physically
meaningful given §2.4's box-mode identification.

**Flagged caveat:** `L_loop ≈ 535 pH` is in the "hundreds of pH" range that `EPR_C_EXTRACTION.md`
§8 identifies as the one place the missing kinetic-inductance modelling (§3 of that document) could
actually matter, since it affects the screening parameter `β_L`. It does not affect the `c`
convergence conclusion above, which is purely a geometric/mesh result.

---

## 3. Bottom line (baseline geometry, `l1=l2=300µm`)

**`c ≈ 8.7×10⁻⁵` rad**, for the physical (readout-resonator) mode, converged with respect to both
mesh refinement (AMR, independent pipeline) and `L_probe` (3-point fit, 2.5% spread). The
coarse-mesh estimate (2.00×10⁻⁴) was off by more than 2× and should not be used.

Outstanding items against the `EPR_C_EXTRACTION.md` validation checklist:
- §6.2 (junction-mesh-density convergence, formal sweep of `p_probe` vs. junction element count)
  was not completed as a standalone check — the SQUID loop's resolution was verified once
  (~3.1µm median) rather than swept.
- The AMR iteration-3 `STable3D` crash was not root-caused; a third, more-refined AMR iteration
  might still shift `c` further, though the `L_probe`-fit agreement with AMR-2 iteration 2 suggests
  this is unlikely to be large.
- Mode 2's box-mode nature suggests the simulation domain (chip size / PEC wall placement) should
  be revisited if a genuine second resonator harmonic is needed.

---

## 4. Geometry change: SQUID arm doubled (`l1=l2=300µm → 600µm`)

To test whether `c` depends on the SQUID loop's own geometry (as opposed to just mesh/`L_probe`
convergence), the design was changed to double the SQUID arm length (`squid_options.l1`/`l2`:
300µm → 600µm each, i.e. 600µm→1200µm total loop-branch length), together with a corresponding
resonator lead adjustment (`readout_1`'s `lead.end_straight`: 575→980→**1180µm**, and
`launch_point1.pos_x`: -350→-275µm) to keep the resonator routed sensibly around the now-longer
arm. All of this was authored interactively in `build_design_visualize.ipynb` and had to be ported
into `build_design.py` (the file the CLI pipeline actually uses) by hand, in three passes, because
two of the notebook's changes were missed on the first two attempts — a real methodological risk
worth flagging: **when the notebook and the production script can diverge, diff every parameter
explicitly rather than skimming for "the part that changed."**

### 4.1 Debugging trail

1. **First port attempt**: only copied `squid_options.l1`/`l2`. Produced **negative** participation
   ratios for the genuine resonator mode (`p_probe = -1.86×10⁻³`, `p_transmon = -1.17×10⁻²` at
   7.87 GHz) — energy fractions that should be strictly non-negative — plus an extra spurious mode
   that didn't exist in the baseline geometry. Root cause, confirmed by direct polygon plotting
   (`geometry_test2.py`, not a GUI screenshot): **the SQUID arm and the resonator's lead-in trace
   were running parallel and within ~20–40µm of each other for over a millimeter**, since only the
   SQUID had been lengthened while the resonator's own lead geometry still assumed the old, shorter
   arm. This produces uncontrolled near-field coupling between the two inductive ports that the
   `EPR_C_EXTRACTION.md` lumped-port participation formula doesn't cleanly decompose, hence the sign
   anomaly.
2. **User fixed `lead.end_straight` (980→1180µm)** in the notebook to route the resonator further
   from the lengthened arm; only *that* one parameter was re-diffed and ported, missing that
   `launch_point1.pos_x` had also changed. Result: still qualitatively wrong (mode 1's participation
   still negative, plus **two** box-like modes now appearing in the capture window instead of one,
   pushing the real transmon mode to position 4).
3. **Full parameter-by-parameter diff** of every cell in `build_design_visualize.ipynb` against
   `build_design.py` (qubit options, launch points, jogs, lead, chip sizing) finally caught the
   `launch_point1.pos_x` discrepancy too. User confirmed the resulting geometry visually as correct.
   With this fully-corrected geometry, mode 1's participation came back **positive**
   (`p_probe = +6.67×10⁻⁷` at `L_probe=10nH`) and stayed positive and internally consistent across
   the full `L_probe` sweep.

**Lesson**: a maximum of one missed parameter is enough to produce electromagnetically-plausible
looking but physically wrong results (finite, sensible-looking - if negative - participation
values, not an obvious crash) — visual/geometric confirmation via a direct polygon plot, not just a
GUI screenshot at full-chip zoom, was what actually caught the problem, and even then it took two
extra iterations of "diff every parameter" to fully resolve.

### 4.2 Mode identification, revisited

With five converged modes near the `Target=7.0GHz` shift, only mode 1 (7.80 GHz) is a genuine
resonator mode by its field plot; modes 2 (9.99 GHz) and 3 (13.6–13.7 GHz) are both box-mode
artifacts of the (now larger, since the chip auto-sizes around the longer arm) simulation domain,
despite mode 3 being labeled `resonator_mode_3` by the frequency-adjacent participation-ranking
logic in `build_and_eigenmode.py`. Mode 4 (13.8 GHz) is the actual transmon-adjacent mode. This
reinforces `EPR_C_EXTRACTION.md` §6.3's warning to identify modes by field pattern/participation,
not by frequency-ordinal position, especially once the chip footprint grows enough to admit more
box modes into the capture window.

### 4.3 `c` result for the new geometry

Using only mode 1 (the confirmed genuine resonator mode), swept the same three `L_probe` points at
the mesh confirmed clean (1,209,618 elements; SQUID-loop region: 329,273 tetrahedra, 2.25µm median
edge length — consistent with the baseline geometry's resolution):

| L_probe | c_meas (mode 1) |
|---|---|
| 10 nH | 1.995×10⁻⁴ |
| 15 nH | 2.067×10⁻⁴ |
| 20 nH | 2.105×10⁻⁴ |

**Fit: `c_∞ = 2.229×10⁻⁴`, `L_loop = 1170.9 pH`.**

Unlike the baseline geometry, the §6.1 "loop inductance negligible" check **fails** here (5.2%
relative spread, vs. 2.5% before) — doubling the arm length roughly doubled `L_loop` too (baseline:
535–677 pH; here: >1.1 nH), which is now a substantial fraction of even the smallest `L_probe`
tested. This also trips the §8 escalation trigger: at `L_loop` in the low-nH range, the missing
kinetic-inductance modeling (SQDMetal is PEC-only) is no longer a safely-negligible approximation,
so `c_∞ = 2.23×10⁻⁴` should be treated as a first estimate under that caveat, not a fully converged
number the way the baseline's `c ≈ 8.7×10⁻⁵` was.

**Conclusion: `c` is geometry-dependent, as expected** (`8.7×10⁻⁵` at `l1=l2=300µm` vs. `2.23×10⁻⁴`
at `l1=l2=600µm`, roughly 2.6× larger for 2× longer arms) — a longer SQUID loop couples more
strongly to the resonator, consistent with `c ∝ M` (mutual inductance) in the derivation in §1. This
run has not yet been cross-checked with AMR the way the baseline was; doing so, and formally
addressing the kinetic-inductance caveat if `c` for this geometry needs to be trusted to the same
precision as the baseline, are the natural next steps.

---

## 5. Reference

- `EPR_C_EXTRACTION.md` (this directory) — full derivation, port setup, and validation
  requirements.
- Minev et al., *Energy-participation quantization of Josephson circuits*, npj Quantum Information
  **7**, 131 (2021), arXiv:2010.00620.
- SQDMetal paper arXiv:2511.01220, Appendix C.
