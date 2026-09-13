# SQUID-array inductive-coupling model

*Draft — `notebooks/design/flux_optimizing/`. Last updated 2026-09-03, sweep still
running (941 / 2400 points). Numbers marked **[partial]** will be refreshed when
`aggregate_sweep.py` runs on the full grid.*

---

## 1. What this computes

The **geometric mutual inductance** `M` between a readout line and a series
SQUID array, reported as **flux per unit line current** through the array,

```
Φ₀ / mA  ≡  ( Σᵢ Φᵢ  /  I_line )  ×  ( 1 mA / Φ₀ )
```

where `Φᵢ` is the magnetic flux threading loop `i`, `I_line` the current in the
readout line, and `Φ₀ = h/2e = 2.0678×10⁻¹⁵ Wb`. `Σᵢ Φᵢ / I` **is** the mutual
inductance `M` (in henries); the `Φ₀/mA` form is just the natural unit for a
flux-tunable element (the array inductance is periodic in `Φ/Φ₀`).

**This is a design-ranking metric, not the true resonator–array coupling.** See
§6.

---

## 2. Physical model

A magnetostatic finite-element solve (PALACE, via SQDMetal):

- Solve `∇×(ν ∇×A) = 0` in the chip + air box, with a **surface-current
  excitation** injecting `I_line` into the readout line's open end.
- All superconductor is a **perfect electric conductor** (`B⊥ = 0`): the metal
  expels field and carries screening currents. The array metal therefore acts as
  a **flux concentrator** — field excluded from the arms is funnelled through the
  loop apertures (this is why thicker arms raise `M`, §7).
- The Josephson junctions are **not** in the model. Each `JJ_gap` is a real
  vacuum break, so every SQUID loop is an *open* ring. At the design junction
  inductance (`L_J = 250 nH`, or even the nominal 15 nH) the screening a closed
  loop would provide is `≈ L_loop / (L_J/2) ≈ 50 pH / 125 nH ≲ 0.05 %` — utterly
  negligible — so the bare open-loop flux **is** the finite-`L_J` answer to
  within noise. Junction inductance only matters below ~1 nH.
- Flux is read out as a Palace `SurfaceFlux` post-processing integral of `B_z`
  over a surface spanning all 23 loop apertures (`loop_aperture_polys`, §5.4).

---

## 3. Geometry

Built by `build_design.build_design(...)`.

### 3.1 SQUID array (`SquidArray` custom qiskit-metal component)

`N_SQUIDS = 23` identical loops tiled along `+x` and wired in series (a tunable
superinductor). Per loop: `plate1` (left contact), `seg_a`/`seg_a_low` and
`seg_b`/`seg_b_low` (the horizontal arms, split by a `JJ_gap` where the junction
would sit), `seg_c` (right vertical), `seg_d` (right contact). Loop `i`'s right
net overlaps loop `i+1`'s left net via an explicit link bar, so the only
galvanic end-to-end path runs through every junction pair.

- **Loop aperture** `i`:  `(2·arm_semi_length + JJ_gap)` wide  ×  `squid_gap`
  tall, left edge at `plate1`'s right edge.
- **Loop pitch** (loops abut):
  `full_w = arm_width + 2·arm_semi_length + JJ_gap + arm_width + segment_length`.
- A `subtract=True` `ground_gap`-wide pocket is carved around the whole chain
  (also clears ground out of every `JJ_gap`, which would otherwise short it).

### 3.2 Readout probe — hairpin (`RouteAnchors`)

A CPW "hairpin" (there-and-back U) straddling the array, replacing the real
λ/4 resonator for the sweep:

- Runs `route_clearance` µm past the **top** of the array, U-turns at the left
  (elbow) end, returns `route_clearance` µm past the **bottom**.
- `route_clearance` is the gap from the wire centreline to the nearest array
  metal, held equal on the top / bottom / left-elbow / termination sides.
- `route_fillet = 20 µm` (< `route_clearance`).
- Terminations: `OpenToGround` (`launch_point1`, current injected here) and
  `ShortToGround` (`launch_point2`, DC return through the ground plane). A
  **U-clip** at the open end ties the two coplanar grounds around the injection
  point so the DC loop closes.
- CPW: `cpw_width = 20 µm`, `cpw_gap = 12.25 µm`.

### 3.3 Array near-end ground bridge (`ground_near_end = True`)

A short positive-metal `Rectangle` ties the array's **left end** (`plate1` of
loop 0, the end nearest the resonator elbow) to the ground plane, bridging the
`ground_gap` pocket. Tested vs fully-floating ends; the difference is ~1 %, and
the user chose the single-ground configuration for the sweep.

### 3.4 Chip

`DesignPlanar`, silicon substrate. Sized to the geometry bounding box + an
**absolute margin** (`chip_margin_um = 250` for the sweep, not the 2.5× factor —
the hairpin field is a localised dipole that decays as `1/r³`). Substrate
thinned to `chip_size_z_um = 200` (magnetostatically inert — silicon has
`µ_r = 1` — so it only trades mesh volume). Air box: `x`/`y` padding 5 % of chip
size, `z` 1.5× substrate thickness above, PEC floor below.

---

## 4. Parameters

### 4.1 Swept axes (`sweep_flux.py::AXES`) — 4 × 5 × 6 × 5 × 4 = **2400 points**

| axis | µm values | nominal | physical meaning |
|---|---|---|---|
| `route_clearance_um` | 30, 33, 36, 40 | 35 | wire-centreline → nearest array metal |
| `segment_length_um`  | 0, 3, 6, 9, 12 | 6 | `seg_d` length = extra inter-loop lead → loop pitch |
| `squid_gap_um`       | 14, 19.2, 24.4, 29.6, 34.8, 40 | 25 | loop **height** (aperture y-extent) |
| `arm_semi_length_um` | 8, 12, 16, 20, 24 | 15 | `seg_a`/`seg_b` length → loop **half-width** |
| `arm_width_um`       | 3, 4.5, 6, 8 | 4 | arm / plate thickness (all `plate1`/`seg_*` widths) |

`arm_width` starts at 3 µm: 2 µm arms do not mesh-converge and are the
lowest-flux corner anyway.

### 4.2 Fixed geometry (`build_design.py` constants)

| constant | value | note |
|---|---|---|
| `N_SQUIDS` | 23 | loops in the chain |
| `JOSEPHSON_GAP_UM` (`JJ_gap`) | 3.0 | the vacuum break where each junction would sit |
| `LINE_WIDTH_UM` (`seg_d` height) | 9.0 | right-contact height |
| `GROUND_GAP_UM` | 3.0 | ground-plane pocket around the chain |
| `JJ_WIDTH_UM` | 2.0 | (junction geometry — unused, no ports) |
| `L_J_NH` / `CJ_J_FF` | 15 / 1.5 | (junction electrical — unused; design target `L_J = 250 nH`) |
| `ARRAY_POS_X_UM`, `ARRAY_POS_Y_UM` | −1290, 0 | array origin |
| `cpw_width`, `cpw_gap` | 20, 12.25 µm | readout CPW |
| `ROUTE_FILLET_UM` | 20 | hairpin corner radius |

---

## 5. Simulation setup (`sweep_flux.py::_attempt`)

### 5.1 Solver

| setting | value | why |
|---|---|---|
| type | PALACE magnetostatic | |
| `solver_order` | 1 | order 2 is 5× the DoF for no ranking benefit |
| Krylov / preconditioner | **GMRES + AMS** | plain CG *stalls* — the 23 open loops are topological handles whose cohomology field modes AMS cannot precondition (residual plateaus ~2e-2). GMRES grinds through them (~300–1600 iters). |
| `MaxSize` (GMRES restart) | 250 | balances restart-thrashing vs memory |
| `solver_tol` | 1e-3 | verified: identical to 1e-7 within 0.06 % (it is genuine convergence, not a loose cap) |
| `solver_maxits` | 1200 | pathological points that hit this are recorded with `converged = False` and an approximate value — flagged, not dropped (~4–5 % of points) |
| `solns_to_save` | 0 | no paraview VTU — only the scalar flux is needed |

### 5.2 Mesh — **uniform 1.5 µm** over the loop band

`fine_mesh_in_rectangle` (2 µm ribbon → set to 1.5) over the loop-aperture bbox +
8 µm pad, and `fine_mesh_components(['array_1'], min_size = 1.5 µm, taper 6 µm)`.
~1.4–3.9 M tets.

Chosen after a convergence study:

| mesh | Φ₀/mA (nominal) | vs converged | cost |
|---|---|---|---|
| uniform 3 µm | 11.9 | −15 % | ~35–120 s |
| **uniform 1.5 µm** | **13.7** | **−1.7 %** (near-flat offset) | **~150–2000 s** |
| boundary-layer ("bl_b") | 13.95 | ref (<1 %) | 30–65 min — infeasible |

The −1.7 % offset is nearly uniform across geometries, so the *ranking* is clean.
For accurate **absolute** `M`, re-run the top configs at bl_b (a ~2 h job).

### 5.3 Boundary conditions

PEC: all metal + air-box far field. `SurfaceCurrent` excitation on a patch of the
readout centre conductor at `launch_point1`. `SurfaceFlux` post-proc over the
integration surface.

### 5.4 Integration surface — `loop_aperture_polys(...)`

23 rectangles, one per loop aperture, in metres, computed from the array
parameters (no need to instantiate `SquidArray`). Each box is grown outward by
`grow_um = 0.3 µm` (retry 1.0, 2.5) so its edges do **not** coincide with the
array-metal edges — an exact coincidence makes GMSH build degenerate zero-length
edges ("Could not create line"). The 0.3 µm sliver sits over PEC where `B_z ≈ 0`,
so it adds ~no flux and, being size-independent, adds no size-dependent bias.

*SQDMetal patch:* `add_integration_area` only handles a single `Polygon`;
`_create_gmsh_geometry_from_shapely_polygons` is monkey-patched to flatten the
23-part `MultiPolygon` into its constituent polygons (→ one `SurfaceFlux` group
covering all loops → Palace returns the summed flux).

---

## 6. Key approximations & caveats

1. **Uniform probe current ≠ mode current.** The hairpin carries a *uniform* DC
   current; the real λ/4 resonator has a current profile (max at the short, zero
   at the open). The true coupling `g` weights the per-loop `M` by the local mode
   amplitude. This sweep gives the **geometric envelope**; layer the mode profile
   on top when the resonator goes back (position the array at the current
   antinode, near the short).
2. **Bare open-loop flux.** Junctions absent; valid because `L_J ≫ L_loop`.
3. **Mesh: −1.7 % low** (uniform offset). Fine for ranking; re-mesh winners for
   absolute `M`.
4. **~4–5 % of points `converged = False`** — the worst near-null-space
   geometries (big loops). Their value is approximate (residual ~1e-3…1e-2);
   the two 12-core gap-fill passes reconverge some.
5. **Array metal *in* the model matters** — dropping it (loops as vacuum
   apertures) gives ~22 % lower flux. It is a real flux-concentration effect, not
   an artefact.

---

## 7. Result: scaling law  **[partial, 941 pts]**

Global multiplicative fit (log-space R² = 0.995, median error 2.1 %):

```
Φ₀/mA  ≈  7.0×10³  ×  route_clearance⁻²·⁹⁶
                    ×  squid_gap⁺⁰·⁶⁵
                    ×  arm_semi_length⁺⁰·⁷¹
                    ×  arm_width⁺⁰·²⁸
                    ×  (segment_length + 3)⁻⁰·⁰⁴        [lengths in µm]
```

| axis | exponent | reading |
|---|---|---|
| `route_clearance` | **−2.96** | flux ∝ 1/distance³ — the hairpin's two antiparallel passes form a magnetic **dipole**; dipole field decays as `1/r³`. By far the strongest lever (~10× per µm vs any other). |
| `arm_semi_length` | +0.71 | loop half-width. Pure-area scaling → +1; the ~+0.7 is **diminishing returns** (non-uniform field; far loop edge sits farther from the wire). |
| `squid_gap` | +0.65 | loop height — same story. |
| `arm_width` | +0.28 | flux-concentrator effect is real but **saturating** (double the arm → ~1.2× flux). |
| `segment_length` | ≈ 0 (−0.04, **preliminary**, still ~seg-0-only sampled) | inter-loop spacing barely affects per-amp flux linkage. |

**Design takeaway:** closing the wire–array gap is the lever that matters. Loop
size helps but saturates; arm width is a minor tweak. Every axis is monotonic →
the optimum is a grid corner (`clr 30, seg 0, gap 40, semi 24, aw 8` →
**49.5 Φ₀/mA** so far, vs the min 3.5). The fit is accurate enough (2 %) to
interpolate without more sweep points.

---

## 8. Cost model  **[partial]**

```
n_tet    ≈  40k + 92k·arm_semi + 108k·arm_width + 22k·squid_gap + 45k·seg     (R² 0.99, µm)
seconds  ≈  122 + 2.26×10⁵·(n_tet · n_ksp / 10¹²)                             (R² 0.65 — n_ksp erratic)
seconds  ∝  n_tet²                                                           (univariate)
```

`n_ksp` (GMRES iterations) ranges 250–1600 and rises with loop size (worse
near-null-space). Per-point wall time: median ~10 min, tail to ~30 min at 8
cores.

---

## 9. Files & reproduction

| file | role |
|---|---|
| `build_design.py` | `build_design(...)` → the full design; `loop_aperture_polys(...)` → integration surface |
| `sweep_flux.py` | one shard of the sweep (grid, mesh, solver, per-point CSV, resume-by-grid-point) |
| `submit_sweep.sh` / `submit_sweep_hi.sh` | SGE array jobs (8-core/44 G main; 12-core/60 G gap-fill) |
| `submit_aggregate.sh` + `aggregate_sweep.py` | merge shards → `sweep_results.csv`, `sweep_colormaps.png`, top configs, `DONE.txt` |
| `count_done.py` | `<n_done> <n_error>` across shard CSVs |
| `viz_sweep_trends.py` | the trends+fits figure (`python viz_sweep_trends.py <dir-with-shard_*.csv>`) |
| `inductance_sim.ipynb` | single-point interactive version (nominal geometry) |

**Cluster run:** `~/circuit_design/flux_optimizing/`. Chain
`qsub submit_sweep.sh` → `qsub -hold_jid <id> submit_sweep_hi.sh` (×2) →
`qsub -hold_jid <id> submit_aggregate.sh`. Currently job `1938674` + held chain.

---

## 10. Status (2026-09-03 20:20 −03)

941 / 2400 points, 0 errors, 42 `converged = False` (~4.5 %). Φ₀/mA range
3.5 → 49.5. ETA ~2026-09-04 midday. `segment_length` still under-sampled.
