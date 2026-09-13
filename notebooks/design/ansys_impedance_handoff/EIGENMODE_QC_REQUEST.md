# Follow-up request: add a lossy port to the existing HFSS eigenmode run

## Who this is for

The same local ANSYS/HFSS agent as `HANDOFF.md`. This is a small, targeted
follow-up to `testC_eigen.py` (the eigenmode run on `full_design_cs14_cl300`),
not a new investigation.

## Why

`testC_eigenmodes.json` / `testC_full_design_fit.json` report `Qc_status:
"not extracted"`. Looking at why: the driven S-sweep on this geometry never
resolved a resonance at all (`|S11|` flat at 0.054±0.002 across 6.5-8.7GHz —
the coupling through the CoupledLineTee's ~15um ground strip is too weak for
a driven adaptive sweep to find), and separately, the eigenmode run
(`testC_eigen.py`) used `add_ports=False` — a bare PEC-walled box with only
silicon's dielectric loss tangent (1.2e-5) as a loss channel. So the reported
`Q=92231` for the 8.546GHz mode is a **dielectric-limited Q**, not a coupling
Q — there is currently no external/coupling Qc from HFSS to compare against
at all.

We hit the exact same "driven notch too shallow to resolve" wall on the
PALACE side, and switched to the same fix proposed below: extract Qc from a
**lossy eigenmode** instead of a driven S-parameter fit.

## The ask -- Part 1: cs14_cl300 (original request, still pending)

Re-run `testC_eigen.py`'s eigenmode setup, but with `add_ports=True` (or
equivalent) and a **lumped RLC boundary of R=50 Ohm** in place of each of the
two CPW feedline ports (same locations `run_ansys_tests.py`'s driven tests
already use). HFSS's Eigenmode solution type does not support true S-parameter
ports (no excitation in an eigenmode solve), but does support a resistive/
lumped RLC boundary condition at a port location — this is the standard
technique (used e.g. throughout pyEPR-style HFSS workflows) for giving an
eigenmode a real loss channel representing an external 50 Ohm environment.

With that resistive boundary in place, the eigenmode solver should return a
COMPLEX eigenfrequency `f = fr + i*fi` for the ~8.5GHz resonator mode (instead
of the current purely-real-up-to-dielectric-loss one). Report, per mode:

- `fr` (real part, Hz)
- `fi` (imaginary part, Hz) -- or however HFSS's own eigenmode report already
  expresses mode "Q" (it may already combine both real R=50 loss AND the
  dielectric tan_delta automatically once the lossy boundary is added -- report
  whatever HFSS gives natively, we can separate the two loss contributions
  numerically on our end if needed)
- HFSS's own reported "Q" for that mode

We'll compute `Qc = fr / (2*fi)` (or however it decomposes) on our end and
compare directly against our own PALACE eigenmode result at the same design
point (`cs14_cl300`): **Qc ≈ 13,562** (2-port, both R=50 Ohm) to **26,966**
(1-port, second port set near-open ~1GOhm) -- we currently don't fully trust
either number and would like an independent check. The driven-S21-fit
reference (also PALACE, less trusted given the shallow-notch problem) was
**Qc ≈ 15,684**.

## The ask -- Part 2 (HIGHER PRIORITY, added 2026-09-09): two production-grid points where we have NO number at all

Since Part 1 was written, we mapped the full 42-point production grid
(`coupling_space_um` x `coupling_length_um`) with the same PALACE lossy-
eigenmode technique, and separately pinned down its reliability limit
precisely: it agrees with a clean analytic reference to 1-2% up to
`Qc≈82,000`, then completely saturates on a ~4.5e7 numerical noise floor by
`Qc≈85,000` -- a hard cliff, confirmed by bisection, not a gradual bias.
**11 of the 42 grid points landed above that ceiling and currently have NO
trustworthy Qc from either PALACE method** (the driven-S21 fit doesn't reach
this regime either -- 8/11 missed the simulation window entirely, the other
3 are noise-dominated fits we'd already independently flagged and excluded).

This is exactly the situation where an independent tool matters most: HFSS's
eigenmode solver uses different underlying numerics (not necessarily the same
absolute noise floor as PALACE's), so it may be able to resolve a real Qc for
these points where PALACE currently cannot.

Two new geometries, generated the same way as `full_design_cs14_cl300.gds`:

- `geometry/full_design_cs22p5_cl250.gds` (`coupling_space_um=22.5,
  coupling_length_um=250`) -- **primary ask**. PALACE eigenmode gives
  `Qc=95,200`, just past our confirmed ceiling -- the most informative single
  point, since it directly tests whether the true Qc is close to that number
  or something else entirely (the same way our own ceiling saturates to a
  fixed ~4.5e7 regardless of the true value once past it).
- `geometry/full_design_cs25_cl100.gds` (`coupling_space_um=25,
  coupling_length_um=100`) -- **secondary/optional, if time allows**. PALACE
  eigenmode gives `Qc=121,515`, further into the saturated zone -- useful for
  seeing whether the trend continues sensibly.

Same method as Part 1 (`add_ports=True`, R=50 Ohm lumped RLC boundary at each
CPW feedline port, report complex `f = fr + i*fi` and HFSS's own "Q") -- no
new technique needed, just two more geometries through the same recipe. No
need to also attempt a driven S-sweep on these -- given the coupling is even
weaker than `cs14_cl300`'s (which already failed), it's extremely unlikely to
resolve, and isn't needed if the eigenmode route works.

## What NOT to worry about

- No need to resolve the driven-sweep shallow-notch problem -- that's a
  separate, likely unresolvable-by-either-tool issue at this coupling
  strength (both PALACE and HFSS hit it independently, which is itself a
  useful data point: it's a physical/numerical limit of S-parameter fitting
  at this Qc, not a bug in either tool).
- Part 1 (cs14_cl300) and Part 2 (the two new points) can be done in either
  order or independently -- they don't depend on each other.
