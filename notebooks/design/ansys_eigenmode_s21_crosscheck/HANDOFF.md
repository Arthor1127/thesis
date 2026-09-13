# Handoff: ANSYS driven-S-parameter cross-check of PALACE's eigenmode Qc method

## Who this is for

The same local ANSYS/HFSS agent as `ansys_impedance_handoff/`. This is a
**separate, self-contained** request — you do not need that directory's
context, just this one.

## Why this geometry, not the real hanger design

We asked ANSYS to check our real hanger design (`full_design_cs14_cl300`,
in `ansys_impedance_handoff/`) before, and its driven S-sweep could not
resolve a resonance at all — `|S11|` stayed flat at 0.054±0.002 across the
whole band. That's a real, physical result (the design is deliberately
weakly coupled, high-Qc), not a bug, but it means that geometry can never
give us an S-parameter-based cross-check.

This directory uses a **different, deliberately simple geometry** instead:
a bare shorted quarter-wave CPW resonator with ONE lumped port at the open
end (not a hanger — no separate feedline, no coupler). Its Qc is *fully
analytic* — no coupling capacitance to estimate — and we control it directly
via the port resistance R, so we can pick R values guaranteed to give a
resolvable notch, unlike the real design.

## Geometry

`geometry/bare_resonator.gds` (see `export_geometry.py` for the generating
script). CPW cross-section `cpw_width=20um, cpw_gap=11.5um` (same as
everywhere else in this project, on 500um silicon `er=11.45`).

- One end (`OTG1`, at x=-1886.295um): an open-circuit stub in the GDS — this
  is where the **lossy port** goes (see below). It's drawn as a plain gap
  etch in the GDS with no taper; place your port/boundary right at this gap.
- Other end (`STG`, at x=+1886.295um): a real short-to-ground.
- Total length ~3.773mm, chosen so the fundamental sits near 8GHz (like the
  rest of this project's resonators) — no physical significance beyond that.

## The ask: driven S11 sweep (this is a 1-port structure — S11, not S21)

At the open end (`OTG1`), terminate the CPW with a **lumped port of
resistance R** (same technique as our other request in
`ansys_impedance_handoff/EIGENMODE_QC_REQUEST.md`, but here as a proper
driven S-parameter PORT, not just a lossy boundary — this geometry has only
one port so there's no S21, only S11 reflection).

Run a driven frequency sweep and fit the S11 resonance for **three separate
R values** (each a separate simulation), matching what we already tested in
PALACE:

| R (Ohm) | Analytic Qc = (pi/4)*(R/50) |
|---------|------------------------------|
| 2e5     | 3144                         |
| 1e6     | 15720                        |
| 5e6     | 78601                        |

For each, report: resonant frequency `fr`, and `Qc` extracted from the S11
notch (standard reflection-resonance fit — e.g. the notch depth/width, or
however your usual driven-Q extraction works). Please also save the raw
`.s2p`/S11-vs-frequency data, same convention as `ansys_impedance_handoff/`.

## What we're comparing against

The formula `Qc = (pi/4)*(R/Z0)` (Z0=49.961 Ohm) is standard cavity
perturbation theory for a resistively-loaded quarter-wave resonator — not
something we're asking you to validate, just the reference point.

Our own PALACE eigenmode result (complex-eigenfrequency method, same R
values, same geometry) already ran:

| R (Ohm) | PALACE eigenmode Qc (loaded_Q) | Qc / analytic |
|---------|--------------------------------|----------------|
| 2e5     | 3841.1                         | 1.222          |
| 1e6     | 18495.4                        | 1.176          |
| 5e6     | 78181.5                        | 0.995          |

If HFSS's driven-S11-fit Qc lands close to the analytic column (not the
PALACE column), that tells us PALACE's eigenmode method has some R-dependent
bias worth chasing further. If it lands close to PALACE's numbers instead,
that's reassuring cross-validation of the eigenmode method from a completely
independent tool and technique (driven vs. eigenmode, HFSS vs. PALACE).
Either outcome is useful — just report the raw fit, no need to interpret it.
