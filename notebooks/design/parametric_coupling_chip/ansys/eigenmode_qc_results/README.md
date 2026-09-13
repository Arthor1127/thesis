# `eigenmode_qc_results/` — HFSS lossy-eigenmode cross-check, handed back

Answers **`../EIGENMODE_QC_REQUEST.md`** (both parts). Read **`SUMMARY.md`** first
— the short version is that HFSS eigenmode **could not resolve the resonator
Qc** for these hanger points on this 16 GB machine (the resonator↔feedline
coupling is below what the adaptive mesh captures, and a fine enough coupler
mesh seed OOM-crashes). What *is* solid: the two new geometries, the
validated 50 Ω port model, and `fr_HFSS ≈ 8.46–8.51 GHz`.

## Files

```
SUMMARY.md                     the written answer  <- start here
RESUME.md                      how to finish this (bigger-RAM machine, or
                                 driven-S21 + coupler seed) with exact commands
qc_vs_palace.png               HFSS eigenmode Qc (noise floor) vs PALACE

eigen_qc_<geom>_result.json    per geometry: resonator-mode complex eigen-
                                 frequency, Q_total, the (noise-floor) Qc,
                                 every mode found (incl. the damped feedline
                                 modes that show the ports work)
eigen_qc.py                    the eigenmode driver (lumped-RLC 50 Ω ports,
                                 optional coupler mesh seed)
export_geometry.py             regenerates the two new GDS
port_validation.py             |S11| check that the RLC port pair terminates
                                 the feedline (-18 dB)
read_solved.py                 recover eigenvalues from a solved project if a
                                 readback crashes
make_figure.py                 rebuilds qc_vs_palace.png
```

`<geom>` ∈ {`cs14_cl300`, `cs22p5_cl250`, `cs25_cl100`}.

## The two new geometries

`../geometry/full_design_cs22p5_cl250.gds` and
`../geometry/full_design_cs25_cl100.gds`, from `export_geometry.py` →
`../source_reference/readout_resonators.py::build_design(coupling_space_um,
coupling_length_um)` — the same generator as the shipped `cs14_cl300.gds`.
Verified: `python export_geometry.py --check` re-generates cs14_cl300 and
bbox-matches the shipped file to 0.0 µm.

## Running

```
conda activate qiskit_qc
python export_geometry.py          # (re)make the two new GDS
python eigen_qc.py                 # all 3 geometries, headless HFSS
python eigen_qc.py cs22p5_cl250    # one
python port_validation.py          # the -18 dB port check
```

Env knobs: `HFSS_CORES` (8), `HVAC` (vacuum box µm), `DF` (MaxDeltaFreq %),
`MAXPASS`/`MINPASS`, `FMIN`, `RES_NMODES`, `SI_TANDELTA` (1.2e-5 lossy /
0 lossless), `PORT_X` (clip domain — leave unset for launchpad-face ports;
`965` moves them to the feedline but they do **not** connect there),
`SEED_COUPLER` (1 to force a fine coupler mesh — needs > 16 GB RAM),
`SEED_MAXLEN`/`SEED_X0..Y1` (seed window), `READ_PASSHIST` (0; the per-pass
history query can OOM the readback), `AEDT_GRPC` (1 = gRPC transport).

## Method

Geometry build = `../ansys_results/run_ansys_tests.build_from_gds`
(zero-thickness PEC from the GDS, grounded backplane, compact box; the hanging
resonator meander is in the ground net — galvanically shorted to ground
through `res_short`; the feedline is the signal conductor). Silicon
`tan δ = 1.2e-5` (`SI_TANDELTA`).

At each domain x-face (x = ±X = the outer launchpad edge, same as the driven
wave ports) a lumped-RLC sheet across **each** CPW gap, conductor→ground,
`R = 100 Ω` per gap → 50 Ω effective. The sheet spans the launchpad-pad edge
(|y| ≈ 125 µm) to the launchpad ground edge (|y| ≈ 275 µm — the
`pad_gap = 150 µm` slot, **not** the 11.5 µm feedline gap), overlapping 3 µm
onto each conductor. `port_validation.py`: `|S11| ≈ −18 dB` into the RLC end
(≈ 98.5 % absorbed). In the eigenmode solve these ports pull the ~16 GHz
feedline mode from Q ≈ 85 000 to Q ≈ 15 000–18 000.

HFSS Eigenmode then returns `f = fr + i·fi` per mode. For the resonator mode
`fi` is at the solver's numerical-noise floor (~1 Hz), so
`Q_total = fr/(2 fi) ≈ Qi ≈ 91 000` (dielectric only) and the extracted
`Qc` (via `1/Qc = 1/Q_total − 1/Qi`, `Qi = 92 231`) is meaningless —
the coupling was not captured. See `SUMMARY.md` / `RESUME.md`.
