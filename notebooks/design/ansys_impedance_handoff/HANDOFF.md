# Handoff: ANSYS cross-check of a CPW impedance discrepancy

## Who this is for

A separate Claude Code agent, running locally (not on the CNEA cluster), with
access to ANSYS (presumably HFSS, possibly Q3D Extractor). Your job is to
build and run the tests described below, then hand the **raw results** back
so the requesting session can analyze them. You do not need to draw
conclusions about *why* any discrepancy exists — just produce clean,
well-documented data.

## Background (read this first)

We are designing a superconducting CPW (coplanar waveguide) circuit in
qiskit_metal, simulated with SQDMetal + PALACE (an open-source FEM EM
solver) on a remote cluster. The CPW cross-section (`cpw_width=20um`,
`cpw_gap=11.5um`, on 500um-thick silicon, `er=11.45`) was chosen so that two
independent analytic formulas both predict `Z0 ≈ 50 Ohm`:

- This project's own `DesignHelper.compute_cpw_geometric_parameters`: **49.96 Ohm**
- SQDMetal's own built-in `CpwParams.calc_impedance` (a different, independent
  conformal-mapping formula that also accounts for the finite substrate
  thickness / grounded backplane): **50.29 Ohm**

But when we simulate the SAME cross-section in PALACE (a driven S-parameter
simulation, extracting `Zin = 50*(1+S11)/(1-S11)` from the reflection
coefficient at a lumped CPW port), we get **Zin ≈ 40.5 Ohm** (`40.45 - 1.66j`),
a **~19% deviation** from the analytic prediction. We spent a full session
ruling out every domain/mesh explanation we could think of:

- Mesh resolution: refined 4x beyond our normal rule — the answer moved from
  33→38→40.5 Ohm as it refined, then **plateaued exactly at 40.5 Ohm** with no
  further movement under 2 more refinement doublings.
- Vacuum height above the chip: varied 1x-12x the substrate thickness — **zero
  effect** (33.3/33.7/33.1 Ohm, pure noise).
- Lateral ground-plane width: varied 0.1x-3.0x the chip size (a 30x range) —
  **zero effect** (identical to 10+ decimal places).
- Material permittivity in the actual PALACE config: confirmed `11.45`,
  an exact match to what both analytic formulas assume.

We then tried a completely different extraction method — capacitance-per-length
via `PALACE_Capacitance_Simulation` (silicon vs. vacuum, giving `Z0` without
ever touching S-parameters or ports) — and got **47.6 Ohm**, only ~5.5% off
the analytic target, which is well within this project's normal
finite-length-fringing residual for this method (other cross-sections tested
earlier gave +9.4% and +5.6% by the same method).

**So**: two independent methods (analytic formulas, capacitance extraction)
agree near 50 Ohm; one method (driven S-parameter `Zin`) disagrees by ~19%,
and is completely insensitive to every mesh/domain parameter we varied. Our
leading hypothesis is that this is a **methodology artifact specific to how
SQDMetal defines its CPW lumped port** (it builds each port as two separate
"fin" elements with opposite field directions, not one simple lumped element
— visible in the generated PALACE `config.json`, see
`source_reference/check_line_impedance.py` and the `Boundaries.LumpedPort`
structure it produces), not a real physical property of the design.

**What we need from ANSYS**: an independent, third EM solver's answer for the
same geometry, ideally using a port/impedance-extraction method that is known
to be robust (HFSS wave ports compute characteristic impedance directly from
the field solution via `Zpi`/`Zpv`/`Zvi`, which is a different and
well-established methodology from PALACE's lumped-port S-parameter approach).
If ANSYS agrees with the ~50 Ohm analytic prediction, that confirms the
PALACE driven-sim port is the outlier. If ANSYS *also* gives ~40 Ohm, the
discrepancy is more fundamental than we think and needs to be revisited from
scratch.

## What's in this directory

```
ansys_impedance_handoff/
├── HANDOFF.md                          (this file)
├── geometry/
│   ├── line_only.gds                   Test A geometry (see below)
│   ├── with_launchpads.gds             Test B geometry
│   └── full_design_cs14_cl300.gds      Test C geometry (optional, lower priority)
├── source_reference/
│   ├── check_line_impedance.py         Exact PALACE setup for Tests A/B (ports, mesh, materials)
│   ├── readout_resonators.py           Exact PALACE setup for Test C (the full design)
│   ├── cpw_impedance_test_build_design.py     Geometry generator used for the capacitance cross-check
│   └── cpw_impedance_test_run_capacitance.py  The capacitance-based Z0 extraction that DID agree with analytic
└── our_results/
    └── impedance_results_summary.json  All our own numbers, in one place
```

The `.gds` files are the **authoritative geometry** — import them directly
into HFSS rather than re-drawing from the Python source. The Python files are
there so you can see exactly what materials/ports/mesh settings PALACE used,
in case you want to mirror them as closely as possible for a fair comparison
(or deliberately do it "the standard ANSYS way" and see if that changes the
answer — either is useful information, just tell us which you did).

## Geometry details (from `impedance_results_summary.json`)

- CPW: `width=20um`, `gap=11.5um` (center trace width / gap to ground each side)
- Substrate: silicon, `er=11.45`, thickness `500um`, effectively lossless
  (`tan_delta≈1.2e-5` in our sim, safe to treat as 0)
- Metal: modeled as a **zero-thickness PEC sheet** in our sim (SQDMetal's
  default "quasi-3D" approximation, not a real finite-thickness conductor).
  If you use a finite-thickness metal (e.g. 100nm aluminum, a common
  realistic choice), that's fine — just note it, since it's a known
  difference from our own setup, not an error.
- Chip size: `4.8mm x 2mm` for Tests A/B, `4.8mm x 4.0mm` for Test C
  (irrelevant to Z0 in principle since we've shown domain size doesn't
  matter, but included for exact reproducibility)

## Tests, in priority order

### Test A (highest priority): bare CPW line impedance

Geometry: `geometry/line_only.gds` — a single straight CPW trace, ~2mm long,
open ends (no launcher pads, no other components).

1. Import into HFSS.
2. Set up the silicon substrate (as above) and vacuum above (any reasonable
   height — we've shown this doesn't matter, but a few mm is standard practice).
3. Terminate each end of the trace with an HFSS **wave port** sized to capture
   the CPW mode's fringing fields (a common rule of thumb: port height/width
   several times the CPW gap). Use HFSS's own **integration line** across the
   gap for the wave port so it computes `Zpi` (or `Zpv`/`Zvi`) directly.
4. Run a driven-modal solve, frequency band **6-10 GHz** (matching our own
   sweep), single frequency point is fine if HFSS's Z0 is frequency-independent
   here (it should be, since there's no resonance in this geometry), but a
   small sweep (e.g. 5-10 points) is good practice to confirm it's flat.
5. **Report**: HFSS's own computed characteristic impedance for the wave port
   mode (`Zpi`, `Zpv`, and/or `Zvi` — report whichever HFSS gives you, ideally
   all three) at at least one frequency, plus the full `S11`/`S21` vs frequency
   (magnitude and phase, or a Touchstone `.s2p` file) so we can independently
   compute `Zin = 50*(1+S11)/(1-S11)` ourselves and compare that too.

### Test B: bare line + launchpad tapers, still no coupler

Geometry: `geometry/with_launchpads.gds` — same CPW trace, but now with a
real launchpad/wirebond-pad taper at each end (a standard impedance-matching
taper from a wide bond pad down to the 20um/11.5um CPW cross-section).

Same procedure as Test A, but now the wave ports sit at the OUTER edge of the
launchpad (the wide bond-pad end), not on the plain CPW trace. This isolates
how much the taper itself contributes to any mismatch. Our own PALACE sim
found `Zin ≈ 29.7 Ohm` here (down from ~50, i.e. the taper alone already
contributes a large chunk of deviation in our sim) — we'd like to know if
ANSYS agrees that the taper itself is a significant contributor, or if it's
still dominated by whatever Test A finds.

### Test C (optional, only if time permits): the full resonator design

Geometry: `geometry/full_design_cs14_cl300.gds` — the complete design: two
launchpads, a feedline, a capacitively-coupled hanger resonator
(`coupling_space=14um`, `coupling_length=300um`).

This is a 2-port driven S21 simulation. Our own PALACE result for this exact
geometry: resonance at **fr = 8.101 GHz**, fitted **Qc = 14078**, **Qr = 533**
(a Khalil et al. 2012 asymmetric hanger-notch fit — see
`our_results/impedance_results_summary.json`'s companion note, or ask for the
fit model if useful), and baseline `Zin ≈ 20.9 + 5.5j Ohm` away from
resonance.

If you run this, sweep **7.5-8.5 GHz** with enough resolution to resolve a
notch that may be quite narrow (our own sim needed ~1600 points across a
750MHz window with adaptive refinement to resolve it cleanly — a Lorentzian
dip with `Qc~14000` at `fr~8.1GHz` has a linewidth of order `fr/Qc ≈ 0.6MHz`,
so if HFSS's own adaptive/interpolating sweep doesn't naturally resolve
something that narrow, a discrete sweep may need many thousands of points or
a tighter window once you've located the dip approximately). Report the full
`S21` magnitude/phase vs frequency (or `.s2p`), plus your own fitted
`fr`/`Q` if you fit it, plus baseline `Zin` away from the resonance the same
way as Tests A/B.

## What to hand back

Please produce, in a directory you create (e.g. `ansys_results/`):

1. **A short written summary** (a few sentences per test) stating what HFSS's
   own field-based impedance calculation gives for Test A, and how it compares
   to 50 Ohm.
2. **Raw S-parameter data** for each test run, as CSV (`freq_hz, S11_mag_db,
   S11_phase_deg, S21_mag_db, S21_phase_deg`) or Touchstone `.s2p` — whichever
   is easier to export from HFSS.
3. **Mesh/convergence info**: number of mesh cells, whether adaptive mesh
   refinement converged, any convergence plot HFSS gives you — just enough
   for us to judge whether the result itself is trustworthy.
4. Optionally, the `.aedt` project file itself if it's not too large, for full
   reproducibility.
5. Note explicitly: (a) whether you used a zero-thickness or finite-thickness
   metal, and if finite, what thickness; (b) exact port/integration-line
   placement choices, since these matter for wave-port field accuracy; (c) any
   deviations you made from this handoff's instructions and why.

Once this is done, hand the results (this whole output directory, or just the
key numbers/files) back to the user, who will pass them to the original Claude
Code session to interpret and decide next steps.
