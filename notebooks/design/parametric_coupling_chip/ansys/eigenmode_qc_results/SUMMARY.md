# HFSS lossy-eigenmode cross-check of PALACE's Qc — result

Request: `../EIGENMODE_QC_REQUEST.md`. Part 1 = `cs14_cl300`; Part 2 =
`cs22p5_cl250` (primary) and `cs25_cl100` (secondary).

## Bottom line

**HFSS's eigenmode solver cannot resolve the resonator Qc for this hanger
geometry on this machine — it hits the same kind of wall PALACE did, for a
different reason.** For every one of the three points, including
`cs14_cl300` where PALACE *has* a number it half-trusts (Qc ≈ 13 562), the
resonator mode's total Q comes back at ≈ 91 000 — the dielectric limit,
unchanged by the R = 50 ports (and ≈ 4·10⁹, `Im f ≈ 1 Hz`, in a lossless
run). The ports themselves work — they damp the ~16.5 GHz *feedline* modes
from Q ≈ 85 000 to Q ≈ 15 000–18 000 — but essentially **no power reaches
them from the resonator mode as HFSS meshes it.**

| geometry | PALACE eigenmode Qc | HFSS resonator `Q_total` | HFSS coupling Qc | HFSS `fr` | HFSS ~16.5 GHz feedline-mode Q |
|---|---|---|---|---|---|
| cs14_cl300  | 13 562 (2-port) / 26 966 (1-port) | 91 002 | **not resolved** (≈ 7·10⁶, noise) | 8.477 GHz | 14 672 |
| cs22p5_cl250 | 95 200 | 91 061 | **not resolved** (≈ 7·10⁶) | 8.487 GHz | 15 568 |
| cs25_cl100  | 121 515 | 91 405 | **not resolved** (≈ 1·10⁷) | 8.512 GHz | 18 359 |

`Q_total` (dielectric `tan δ = 1.2e-5`, so the dielectric-only `Qi ≈ 92 231`
from the existing no-port run) comes back at ≈ 91 000 for **all three** — i.e. the R = 50 ports add
essentially no loss to the resonator mode, so the extracted coupling Qc
(`1/Qc = 1/Q_total − 1/Qi`) is a ratio of two nearly-equal numbers and is
meaningless (7·10⁶ ± a lot). A lossless (`SI_TANDELTA=0`) run gives the same
story directly: `Im f ≈ 1 Hz` → Qc ≈ 4·10⁹. Full numbers in the
`eigen_qc_*_result.json`.

The **feedline** modes near 16.5 GHz, by contrast, are cleanly damped from
Q ≈ 85 000 (no ports) to Q ≈ 15 000–18 000 — the ports are doing their job;
it is specifically the resonator→feedline coupling that the mesh misses.

## Why — and why it is still a useful data point

The resonator couples to the feedline through the CoupledLineTee's ~15 µm
gap (metal-to-metal ≈ 37 µm for cs14). The eigenmode adaptive mesh refines
where the mode's field is strong — inside the λ/4 resonator — and leaves that
gap coarse, so the eigensolver sees an essentially isolated resonator and its
`Im f` collapses to numerical noise. This is the **same failure mode the
driven Test C sweep hit** (`../ansys_results/testC_full_design_fit.json`:
"the driven adaptive mesh converges on the smooth feedline transmission and
never resolves the weak coupling capacitance across the CoupledLineTee's
~15 µm ground strip"). Now the *eigenmode* adaptive does the same thing.

Forcing a fine mesh in that gap (a floating z = 0 sheet + length-mesh
operation over the coupler window) **does** change the picture — but a mesh
fine enough to resolve the coupling (~5 µm target, ~10⁵ tets just in the
seed, growing to > 2.5·10⁵ over the adaptive passes) **exhausts the 16 GB
RAM and the eigensolver OOM-crashes**. A coarser seed that fits (~8 µm)
still leaves the resonator mode at Qc ≈ 4·10⁹. There is no seed density that
both fits in 16 GB and captures the coupling on this machine.

So: an independent tool, run two different ways (this eigenmode attempt and
the earlier driven S-sweep), **independently confirms that the coupling at
these design points is below what an S-parameter / eigenmode extraction can
resolve without a very fine, expensive coupler mesh** — which is consistent
with PALACE's own ~4.5·10⁷ saturation floor being a numerical artifact, not a
physical Qc. It does *not* give you the "close to 95 200 or something else
entirely?" answer you were hoping for; getting that needs either a
bigger-RAM machine or the driven-S21-with-coupler-seed route (details in
`RESUME.md`).

## What IS solid

- **The two new geometries.** `../geometry/full_design_cs22p5_cl250.gds` and
  `../geometry/full_design_cs25_cl100.gds`, from `export_geometry.py` →
  `../source_reference/readout_resonators.py::build_design(cs, cl)` — the
  same generator as the shipped `cs14_cl300.gds`. Verified:
  `python export_geometry.py --check` re-generates cs14_cl300 and
  bbox-matches the shipped file to 0.0 µm.

- **The lumped-RLC 50 Ω port model.** `eigen_qc.add_feedline_rlc`: two RLC
  sheets per port face (one per CPW gap, 100 Ω each → 50 Ω conductor-to-
  ground). `port_validation.py` (wave port one end, RLC pair the other, on
  `with_launchpads.gds`): **`|S11| ≈ −18 dB`** — ≈ 98.5 % absorbed, a
  well-matched 50 Ω termination. In the eigenmode solve, with the ports on
  the outer launchpad faces, they pull the ~16 GHz feedline mode from
  Q ≈ 85 000 (no ports) to Q ≈ 17 000. The ports do their job; it is the
  resonator→feedline coupling that the mesh misses.
  (Note: a variant with the ports moved inboard to the clean 11.5 µm feedline
  CPW, `PORT_X=965`, did **not** connect — all modes stayed dielectric-
  limited — so the launchpad-face ports are the ones to use.)

- **`fr_HFSS` = 8.477 / 8.487 / 8.512 GHz** (cs14 / cs22p5 / cs25),
  consistent with the earlier no-port eigenmode run (8.55 GHz) and ~4–5 %
  above PALACE's 8.10 GHz (driven) / 8.25–8.28 GHz (eigenmode). The offset is
  the closed-PEC-box artifact that
  `../source_reference/readout_resonators.py` already documents for this
  no-real-environment setup.

## Files

- `eigen_qc_<geom>_result.json` — per geometry: the resonator mode's complex
  eigenfrequency, `Q_total_hfss`, the (noise-floor) `Qc`, every mode found.
- `eigen_qc.py` — the eigenmode driver (lumped-RLC ports; `SEED_COUPLER`,
  `PORT_X`, `HVAC`, `DF`, `MAXPASS` … env knobs; see `RESUME.md`).
- `export_geometry.py` — regenerates the two new GDS.
- `port_validation.py` — the `|S11|` port-absorption check.
- `read_solved.py` — recover eigenvalues from a solved project if a readback
  crashes.
- `make_figure.py` — `qc_vs_palace.png`.
- `RESUME.md` — the two fallback routes (bigger RAM; driven S21 + coupler
  seed) with exact commands.
