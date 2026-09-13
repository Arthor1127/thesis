# Direct B-field cross-check of the port-EPR `c` — single-SQUID, 600 µm arm

**Question:** does the direct magnetic-field method developed for the SQUID
array (`../squid_array_coupling/extract_M_bfield.py`, `C_CALCULATION.md` §2b)
reproduce the established port-EPR `c` for the single-SQUID design, where a
trusted reference exists?

**Answer: yes, to ~5%.** `c_bfield ≈ 2.0×10⁻⁴ rad` vs the port-EPR
`c ≈ 2.0–2.1×10⁻⁴` (`C_EXTRACTION_RESULTS.md` §4.3), L_probe-independent, and
`M²/(L_probe L_r) ≈ p_probe`. This is the one clean validation of the B-field
method against a ground truth.

---

## 1. Setup

- **Geometry:** the current `build_design.py` — `CircTransmonSQUID` with
  `squid_options.l1 = l2 = 600µm`, `end_straight = 1180µm`,
  `launch_point1.pos_x = -275µm` (the "SQUID arm doubled" geometry of §4). The
  readout resonator is **identical to the SQUID array's**: a hairpin fold with
  legs at `y = ±50µm`, `x` from -1330 to -150µm (verified by polygon dump —
  same launch points, jogs, `cpw_width=20µm`, `cpw_gap=12.25µm`,
  `res_length=3.74575mm`).
- **Re-run:** the old `lprobe_sweep_l1_600um_v2` runs did not keep ParaView
  fields, so the eigenmode solve was re-run at `L_probe = 10 / 15 / 20 nH`
  with fields saved (`build_and_eigenmode.py`, `--number-of-freqs 3
  --starting-freq 7.0e9`), on the cluster in
  `~/circuit_design/parametric_coupling_test/bfield_600um/` (jobs 1938553–9).
  Resonator mode = m=1 @ 7.8019 GHz. `p_probe` reproduces the sweep exactly
  (4.77e-7 at 15 nH → port-EPR `c = 2.067e-4`).

## 2. The SQUID loop for the fluxoid — NOT the full arm span

From `circle_transmon_squid.py`: each branch runs pad chord → feed of length
`l1` → **junction** → tail of length `l2` → **into the ground plane** (the
tails do not rejoin). The loop that both junctions bound, and whose flux
enters `φ_a − φ_b`, is `loop_opening`: **pad chord (x ≈ −89µm) to the
junctions (x ≈ −689µm), `|y| < 8.7µm`**. The source comment calls it "the
actual SQUID loop area." The region beyond the junctions (`upper_void`,
x < −689µm) is a ground-plane void — its flux does **not** link the circuit.

**Integrating the full arm span (including `upper_void`) over-counts the flux
~1.8×** — the void sits nearer the resonator's high-current turn and carries
*more* flux than `loop_opening`, but none of it is the SQUID loop flux. Caught
by `port-I.csv`: `Φ_link = (L_probe + L_loop) I_probe ≈ 3.9×10⁻¹³ Wb`,
L_probe-independent; the full-span integral gave ~7×10⁻¹³.

## 3. Flux-surface height

`Φ(z)` over `loop_opening` from a pyvista z-slice is not flat:

| z (µm) | Φ (×10⁻¹³ Wb) | note |
|---|---|---|
| 0.5 | 3.11 | near-surface: cells straddle the zero-thickness PEC sheets, `B_z` diluted by the `B=0` on-metal value → under-reads |
| 1.0 | 3.29 | |
| 1.5–1.75 | 3.45–3.52 | loop flux resolved |
| **2.0** | **3.90** | ≈ `Φ_link` from `port-I.csv` |
| 2.5 | 5.03 | now catching resonator-leg flux curving through above the loop → over-reads |
| 3.0 | 6.01 | |
| 5.0 | 8.80 | |

The physical loop flux sits at the crossover, `z ≈ 1.5–2µm`,
`Φ ≈ 3.5–3.9×10⁻¹³ Wb`. (For the array's wider ~25µm slots this z-systematic
is milder.)

## 4. Result

`L_r = 2 U_m_tot / I_r² = 0.94 nH` (matches a hand estimate for the λ/2 CPW
antinode; same resonator → same as the array's 0.87 nH). Units all pass:
`∫ U_m dV = domain-E E_mag` exactly, `U_m` density `= (|Br|²+|Bi|²)/(2µ0)`.

| L_probe | c_portEPR (ref) | c_bfield (z=2µm) | c_bfield / c_portEPR | M_bfield | M from √(p L_probe L_r) |
|---|---|---|---|---|---|
| 10 nH | 1.995e-4 | 1.96e-4 | 0.98 | 2.47 pH | 2.51 pH |
| 15 nH | 2.067e-4 | 1.99e-4 | 0.96 | 2.50 pH | 2.56 pH |
| 20 nH | 2.105e-4 | 2.00e-4 | 0.95 | 2.52 pH | 2.57 pH |

- **c_bfield ≈ c_portEPR to ~5%.**
- **Both L_probe-independent** (as `c ∝ M`, geometric, must be).
- **`M²/(L_probe L_r) ≈ p_probe`** — ratio 0.90–0.96 at z=2µm — the two methods
  are internally consistent, not just numerically close.
- Reading the flux at z ≈ 1µm instead gives c_bfield ≈ 1.7e-4 (0.8× ref); the
  flux-surface choice is the ~±20% systematic of the B-field method for a
  narrow loop.

**Bottom line:** the direct B-field method reproduces the port-EPR
`c ≈ 2×10⁻⁴ rad` for the single-SQUID 600µm design. The port-EPR method (which
this directory established) and the B-field method (developed for the array)
agree where a reference exists.

---

## 5. Files

| file | role |
|---|---|
| `extract_M_bfield.py` | the B-field method, adapted to the single SQUID (`loop_opening` region, z-scan) |
| `bfield_data/single_600um_Lprobe{10,15,20}_z2.npz` | raw arrays (Φ, I_cpw, M, L_r, c, p_probe, …) |
| `build_and_eigenmode.py` / `build_design.py` | the 600µm geometry + eigen solve |
| `C_EXTRACTION_RESULTS.md` §4 | the port-EPR reference (`c_∞ = 2.23e-4`, 3-point L_probe fit) |

Cluster: `~/circuit_design/parametric_coupling_test/bfield_600um/` (jobs
1938553 solve, 1938558/9 extract).
