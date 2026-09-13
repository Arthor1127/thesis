# Plot data — JJ-array eigenmode study

All arrays are plain NumPy; load with `np.load(path, allow_pickle=True)`.
Everything here is regenerable: `python export_plot_data.py` (dispersion),
`export_field_data.py` on the cluster (fields), and the geometry snippet in the
session notes.

Array lengths studied: **N = 10, 23, 40, 100** junctions (N+1 islands).
Fixed input: **L_J = 1.2 nH** per junction, C_J = 2 fF lumped on each port,
50 µm ground-plane pocket. Palace eigenmode solve, no EPR.

---

## `dispersion_data.npz` — panels of `dispersion_fit.png`

Fit model (open / free-free boundary conditions):

    omega(k)^2 = (2/L_J) (1 - cos k) / (C_g + 2 C_Jeff (1 - cos k))
    k_m = (m + delta) * pi / (N + 1),   m = 1 .. N

**Fit results (scalars):**

| key | value | meaning |
|---|---|---|
| `L_J_H` | 1.2e-9 | per-junction inductance (fixed) |
| `C_g_F` | ~9.46e-16 | island → ground capacitance |
| `C_Jeff_F` | ~9.38e-15 | effective per-junction capacitance (2 fF port + geometry) |
| `delta` | ~0.021 | fitted mode-index offset (≈0 ⇒ ideal free-free) |
| `rms_rel` | ~0.030 | rms relative residual of the global fit |
| `Z_array_ohm` | ~1126 | sqrt(L_J / C_g) |
| `f_plasma_GHz` | ~47.4 | 1/(2π sqrt(L_J C_Jeff)) |
| `f_bandedge_GHz` | ~46.8 | dispersion at k = π |
| `NS` | [10 23 40 100] | |

**Per N (prefix `n{N}_`):**

| key | shape | meaning |
|---|---|---|
| `n{N}_m` | (M,) | mode index 1..M (array modes only) |
| `n{N}_f_GHz` | (M,) | eigenfrequency |
| `n{N}_k_over_pi` | (M,) | (m + delta)/(N+1) — x-axis of the dispersion panel |
| `n{N}_p_abs_sum` | (M,) | Σ_j |p_j| (≈1 for array modes) |
| `n{N}_p_per_junction` | (M, N) | per-junction participation (Palace port-EPR) |
| `n{N}_resid_pct` | (M,) | 100 (f_sim − f_fit)/f_fit |
| `n{N}_box_f_GHz` | (B,) | dropped box/package modes (Σ|p_j|<0.5); only N=10 has any (4, at 78–115 GHz) |

**Fit curve:** `fit_k_over_pi` (600,), `fit_f_GHz` (600,) — the smooth line.

Also `dispersion_data.csv`: flat table `N, m, f_GHz, k_over_pi, p_abs_sum, resid_pct` (51 rows).

---

## `field_data_n{N}.npz` — the per-mode E-field maps (`mode_fields/n{N}/*.png`)

FEM solution sampled onto a fixed grid in a horizontal plane `z_slice_um` above
the metal, over the array region, plus a dense chain-axis cut.

| key | shape | meaning |
|---|---|---|
| `x_um`, `y_um` | (nx,), (ny,) | grid axes (µm), nx≈600, ny≈160 |
| `z_slice_um` | scalar | height of the plane above the metal (2 µm) |
| `f_GHz` | (M,) | mode frequencies (same order as the arrays below) |
| `Emag` | (M, ny, nx) | \|E\| = sqrt(\|E_re\|²+\|E_im\|²) in V/m |
| `Ez_re`, `Ex_re` | (M, ny, nx) | real part of the vertical / axial E component |
| `mask` | (ny, nx) bool | grid points inside the FEM domain |
| `prof_x_um` | (nx,) | chain-axis cut coordinate |
| `prof_Ez_re` | (M, nx) | raw real E_z along y=0 at z_slice — has per-island near-field ripple |
| `island_x_um` | (N+1,) | island centre positions along the chain |
| `island_Ez_re` | (M, N+1) | real E_z averaged over each island footprint ∝ island charge — **the clean standing-wave pattern** (mode m ⇒ m sign changes) |
| `array_x_um` | (2,) | metal extent of the chain along x |
| `pitch_um` | scalar | island pitch (pad_length + pad_gap) |

For a standing-wave / mode-shape figure use `island_x_um` vs `island_Ez_re[i]`;
for a 2D field map use `Emag[i]` (or `Ez_re[i]`) with `extent=[x_um[0], x_um[-1],
y_um[0], y_um[-1]]`. Modes with `f_GHz` above the array band (only N=10: the last
2 of 12, ~78/80 GHz) are box modes.

Mode `i` (0-based) ↔ `mode{i+1}.png` ↔ frequency `f_GHz[i]`.

---

## `geometry_data.npz` — the layout sketch (`josephson_array_preview.png`)

Units: **mm**. Per N (prefix `n{N}_`):

| key | shape | meaning |
|---|---|---|
| `n{N}_islands` | (N+1, 5, 2) | island rectangles, polygon vertices (x,y) |
| `n{N}_ground_pocket` | (1, K, 2) | the subtract=True ground-plane pocket polygon |
| `n{N}_junction_lines` | (N, 2, 2) | junction bridge segments, endpoints (x,y) |
| `n{N}_jj_width_um` | scalar | transverse width assigned to each junction |
