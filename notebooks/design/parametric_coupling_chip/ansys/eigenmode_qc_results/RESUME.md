# RESUME — how to finish the EIGENMODE_QC_REQUEST cross-check

Answering `../EIGENMODE_QC_REQUEST.md` (Part 1 cs14_cl300 + Part 2
cs22p5_cl250 / cs25_cl100). **`SUMMARY.md` is the delivered answer** — HFSS
eigenmode could not resolve the coupling Qc on this 16 GB machine. This file
is the two routes to actually get the number, for whoever picks it up.

Also-known operational hazard: after a run of large eigenmode solves that
OOM-crash, the AEDT COM/RPC layer wedges (`_oproject.InsertDesign(...)` →
`com_error (-2147023170, 'The remote procedure call failed.')`) and only a
machine **reboot** clears it (`-RegServer`, killing processes, clearing
`*.lock` do not). Pace the batch and keep peak RAM well under the limit.

## Done and solid (no HFSS needed)

- **New geometries generated & verified.** `../geometry/full_design_cs22p5_cl250.gds`
  and `../geometry/full_design_cs25_cl100.gds`, via `export_geometry.py` →
  `../source_reference/readout_resonators.py::build_design(cs, cl)`.
  `python export_geometry.py --check` re-generates cs14_cl300 and bbox-matches
  the shipped file to 0.0 µm — the generator is correct.
- **Lumped-RLC 50 Ω port model built and validated.** `eigen_qc.add_feedline_rlc`:
  two RLC sheets per port face (one per CPW gap, 100 Ω each → 50 Ω effective
  conductor-to-ground), grown 3 µm onto each conductor. `port_validation.py`
  on `with_launchpads.gds` (wave port one end, RLC pair the other):
  **`|S11| ≈ −18 dB`** — ≈ 98.5 % absorbed, a well-matched 50 Ω termination.
  A bug was found and fixed here: `_gap_edges_at_face` used to fall back to
  the 11.5 µm feedline gap at the launchpad face (real launchpad gap 150 µm),
  so the sheet floated and did nothing (`|S11| = 0 dB`).

## The physics obstacle (final: ~12 eigenmode runs, incl. 3 clean ones post-reboot)

HFSS Eigenmode with the R = 50 launchpad-face ports **damps the feedline
modes correctly** — the ~16.5 GHz feedline mode goes from Q ≈ 85 000
(no ports) to Q ≈ 15 000–18 000 for all three geometries. But **the
resonator mode's coupling to the ports is not captured**: its `Q_total`
comes back at ≈ 91 000 (= the dielectric `Qi`) for all three, and in a
lossless run `Im f ≈ 1 Hz` → Qc ≈ 4·10⁹ (eigensolver noise floor). The
adaptive mesh refines inside the resonator (strong field) and leaves the
~15–37 µm CoupledLineTee gap coarse — the same failure the driven Test C
sweep hit.

A coupler mesh seed (floating z = 0 sheet + `assign_length_mesh` over the
coupler window) fine enough to resolve that gap (`SEED_MAXLEN` ≈ 5 µm) either
OOM-crashes or, at a coarser `SEED_MAXLEN` ≈ 8 µm that fits, still returns
Qc ≈ 4·10⁹. No seed density both fits in 16 GB and captures the coupling.

`fr_HFSS` = 8.477 / 8.487 / 8.512 GHz (cs14 / cs22p5 / cs25) — consistent,
~4–5 % above PALACE (8.10 driven / 8.25–8.28 eigenmode), the closed-PEC-box
artifact.

## To actually get the coupling Qc — two routes

### Route 1 — coupler-seeded eigenmode on a bigger-RAM machine (≥ 24 GB)

```
conda activate qiskit_qc
cd circuit_design/ansys_impedance_handoff/eigenmode_qc_results
export HFSS_CORES=8 HVAC=250 DF=0.2 MAXPASS=16 MINPASS=8 RES_NMODES=3 \
       FMIN=6.0 SEED_COUPLER=1 SEED_MAXLEN=5.0 READ_PASSHIST=1 SI_TANDELTA=0
python -u eigen_qc.py cs14_cl300        # validate: PALACE has ~13.5k here
#   Qc ~ 1e4-3e4 near 8.4 GHz and stable across passes  -> method works,
#   run cs22p5_cl250 + cs25_cl100, rerun make_figure.py, update SUMMARY.
python read_solved.py <geom> 3          # if a post-solve readback still crashes
```
`AEDT_GRPC=1` forces gRPC transport (COM fallback). `PORT_X` clips the domain
but **do not** use `PORT_X=965` — the ports don't connect on the bare
feedline gap there; leave it unset for launchpad-face ports.

### Route 2 — driven S21 + coupler seed + fine sweep (works on 16 GB, heavier)

Reuse `../ansys_results/run_ansys_tests.py` Test C with `TESTC_SEED_BOX` set
to the coupler window (`-180,130,-75,15,25,5` style) and a fine **discrete**
sweep in a ±few-linewidth band centred on `fr ≈ 8.4 GHz` (the bareres trick:
known `fr` → narrow window → fine mesh). An S21 notch fit is more robust than
a fragile eigenmode `Im f`.

## Files (delivered)

- `SUMMARY.md` — the answer. `README.md` — method / how-to-run.
- `eigen_qc_<geom>_result.json` — the 3 clean runs (resonator `Q_total`
  ≈ 91 000, noise-floor `Qc`, the damped feedline modes, `fr`).
- `eigen_qc_cs14_cl300_SEEDED_noisefloor.json` — the coarse-seed run, kept to
  show the seed did not help either.
- `eigen_qc.py` — the driver. `export_geometry.py` — the two new GDS.
- `port_validation.py` — the −18 dB port check. `read_solved.py` — readback
  recovery. `make_figure.py` — `qc_vs_palace.png`. `_all_results.json` — the
  3 runs aggregated.
