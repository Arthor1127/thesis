# Interdigital transmon-coupler geometry gallery

Regenerate with `conda run -n qcg-quantum-design python render_gallery.py`.

Every panel: **left** = full two-transmon design, **right** = zoom on the Q1 junction. Blue = metal, red hatch = ground-plane cutout. A `[!]` in a title flags a variant whose quick geometry check is off (want: 2 metal nets, 1 connected ground piece, comb-to-comb gap > `finger_gap`).

**Baseline** ([00_baseline.png](00_baseline.png)): pad_radius 90 um, pad_gap 10 um, finger_count 8, finger_length 80 um, finger_gap 4 um, spine_width 12 um, lead_length 40 um, lead_blend_radius 40 um.

Each sweep below also has a stacked contact sheet (`<sweep>.png`).

## Qubit spacing (lead_length: spine-to-pad gap)

`lead_length_um` swept over [14, 24, 40, 60, 90]. Contact sheet: [01_qubit_spacing.png](01_qubit_spacing.png)

- `lead_length = 14 um` -> spacing 310 um, nets 2, min gap 4.0 um, ground 1  ([lead_length_um_14.png](01_qubit_spacing/lead_length_um_14.png))
- `lead_length = 24 um` -> spacing 330 um, nets 2, min gap 4.0 um, ground 1  ([lead_length_um_24.png](01_qubit_spacing/lead_length_um_24.png))
- `lead_length = 40 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([lead_length_um_40.png](01_qubit_spacing/lead_length_um_40.png))
- `lead_length = 60 um` -> spacing 402 um, nets 2, min gap 4.0 um, ground 1  ([lead_length_um_60.png](01_qubit_spacing/lead_length_um_60.png))
- `lead_length = 90 um` -> spacing 462 um, nets 2, min gap 4.0 um, ground 1  ([lead_length_um_90.png](01_qubit_spacing/lead_length_um_90.png))

## Tangent-flare radius where the neck leaves the pad

`lead_blend_radius_um` swept over [10, 20, 40, 60]. Contact sheet: [02_tangent_flare_radius.png](02_tangent_flare_radius.png)

- `lead_blend_radius = 10 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([lead_blend_radius_um_10.png](02_tangent_flare_radius/lead_blend_radius_um_10.png))
- `lead_blend_radius = 20 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([lead_blend_radius_um_20.png](02_tangent_flare_radius/lead_blend_radius_um_20.png))
- `lead_blend_radius = 40 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([lead_blend_radius_um_40.png](02_tangent_flare_radius/lead_blend_radius_um_40.png))
- `lead_blend_radius = 60 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([lead_blend_radius_um_60.png](02_tangent_flare_radius/lead_blend_radius_um_60.png))

## Neck width at the pad edge

`lead_width_um` swept over [6, 12, 22, 36]. Contact sheet: [03_neck_width.png](03_neck_width.png)

- `lead_width = 6 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([lead_width_um_6.png](03_neck_width/lead_width_um_6.png))
- `lead_width = 12 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([lead_width_um_12.png](03_neck_width/lead_width_um_12.png))
- `lead_width = 22 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([lead_width_um_22.png](03_neck_width/lead_width_um_22.png))
- `lead_width = 36 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([lead_width_um_36.png](03_neck_width/lead_width_um_36.png))

## Ground-moat rectangle corner radius

`etch_corner_radius_um` swept over [3, 12, 24, 44]. Contact sheet: [04_moat_corner_radius.png](04_moat_corner_radius.png)

- `etch_corner_radius = 3 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([etch_corner_radius_um_3.png](04_moat_corner_radius/etch_corner_radius_um_3.png))
- `etch_corner_radius = 12 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([etch_corner_radius_um_12.png](04_moat_corner_radius/etch_corner_radius_um_12.png))
- `etch_corner_radius = 24 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([etch_corner_radius_um_24.png](04_moat_corner_radius/etch_corner_radius_um_24.png))
- `etch_corner_radius = 44 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([etch_corner_radius_um_44.png](04_moat_corner_radius/etch_corner_radius_um_44.png))

## Fillet where the moat meets each transmon pocket

`etch_blend_radius_um` swept over [2, 10, 20, 35]. Contact sheet: [05_moat_pocket_cusp_fillet.png](05_moat_pocket_cusp_fillet.png)

- `etch_blend_radius = 2 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([etch_blend_radius_um_2.png](05_moat_pocket_cusp_fillet/etch_blend_radius_um_2.png))
- `etch_blend_radius = 10 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([etch_blend_radius_um_10.png](05_moat_pocket_cusp_fillet/etch_blend_radius_um_10.png))
- `etch_blend_radius = 20 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([etch_blend_radius_um_20.png](05_moat_pocket_cusp_fillet/etch_blend_radius_um_20.png))
- `etch_blend_radius = 35 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([etch_blend_radius_um_35.png](05_moat_pocket_cusp_fillet/etch_blend_radius_um_35.png))

## Ground-moat width (dielectric gap to ground)

`cap_gap_ground_um` swept over [3, 6, 12, 20]. Contact sheet: [06_moat_width.png](06_moat_width.png)

- `cap_gap_ground = 3 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([cap_gap_ground_um_3.png](06_moat_width/cap_gap_ground_um_3.png))
- `cap_gap_ground = 6 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([cap_gap_ground_um_6.png](06_moat_width/cap_gap_ground_um_6.png))
- `cap_gap_ground = 12 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([cap_gap_ground_um_12.png](06_moat_width/cap_gap_ground_um_12.png))
- `cap_gap_ground = 20 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([cap_gap_ground_um_20.png](06_moat_width/cap_gap_ground_um_20.png))

## Finger-root / spine-corner fillet radius

`fillet_radius_um` swept over [0, 3, 6, 10]. Contact sheet: [07_finger_root_fillet.png](07_finger_root_fillet.png)

- `fillet_radius = 0 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([fillet_radius_um_0.png](07_finger_root_fillet/fillet_radius_um_0.png))
- `fillet_radius = 3 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([fillet_radius_um_3.png](07_finger_root_fillet/fillet_radius_um_3.png))
- `fillet_radius = 6 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([fillet_radius_um_6.png](07_finger_root_fillet/fillet_radius_um_6.png))
- `fillet_radius = 10 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([fillet_radius_um_10.png](07_finger_root_fillet/fillet_radius_um_10.png))

## Number of interdigital fingers (coupling strength)

`finger_count` swept over [4, 6, 8, 12, 16]. Contact sheet: [08_finger_count.png](08_finger_count.png)

- `finger_count = 4` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([finger_count_4.png](08_finger_count/finger_count_4.png))
- `finger_count = 6` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([finger_count_6.png](08_finger_count/finger_count_6.png))
- `finger_count = 8` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([finger_count_8.png](08_finger_count/finger_count_8.png))
- `finger_count = 12` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([finger_count_12.png](08_finger_count/finger_count_12.png))
- `finger_count = 16` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([finger_count_16.png](08_finger_count/finger_count_16.png))

## Finger length / overlap (coupling strength)

`finger_length_um` swept over [40, 60, 80, 110, 150]. Contact sheet: [09_finger_length.png](09_finger_length.png)

- `finger_length = 40 um` -> spacing 322 um, nets 2, min gap 4.0 um, ground 1  ([finger_length_um_40.png](09_finger_length/finger_length_um_40.png))
- `finger_length = 60 um` -> spacing 342 um, nets 2, min gap 4.0 um, ground 1  ([finger_length_um_60.png](09_finger_length/finger_length_um_60.png))
- `finger_length = 80 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([finger_length_um_80.png](09_finger_length/finger_length_um_80.png))
- `finger_length = 110 um` -> spacing 392 um, nets 2, min gap 4.0 um, ground 1  ([finger_length_um_110.png](09_finger_length/finger_length_um_110.png))
- `finger_length = 150 um` -> spacing 432 um, nets 2, min gap 4.0 um, ground 1  ([finger_length_um_150.png](09_finger_length/finger_length_um_150.png))

## Gap between interleaved fingers (coupling strength)

`finger_gap_um` swept over [3, 4, 6, 10]. Contact sheet: [10_finger_gap.png](10_finger_gap.png)

- `finger_gap = 3 um` -> spacing 362 um, nets 2, min gap 3.0 um, ground 1  ([finger_gap_um_3.png](10_finger_gap/finger_gap_um_3.png))
- `finger_gap = 4 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([finger_gap_um_4.png](10_finger_gap/finger_gap_um_4.png))
- `finger_gap = 6 um` -> spacing 362 um, nets 2, min gap 6.0 um, ground 1  ([finger_gap_um_6.png](10_finger_gap/finger_gap_um_6.png))
- `finger_gap = 10 um` -> spacing 362 um, nets 2, min gap 6.0 um, ground 1  ([finger_gap_um_10.png](10_finger_gap/finger_gap_um_10.png))

## Transmon pad radius

`pad_radius_um` swept over [60, 90, 120]. Contact sheet: [11_pad_radius.png](11_pad_radius.png)

- `pad_radius = 60 um` -> spacing 302 um, nets 2, min gap 4.0 um, ground 1  ([pad_radius_um_60.png](11_pad_radius/pad_radius_um_60.png))
- `pad_radius = 90 um` -> spacing 362 um, nets 2, min gap 4.0 um, ground 1  ([pad_radius_um_90.png](11_pad_radius/pad_radius_um_90.png))
- `pad_radius = 120 um` -> spacing 422 um, nets 2, min gap 4.0 um, ground 1  ([pad_radius_um_120.png](11_pad_radius/pad_radius_um_120.png))

