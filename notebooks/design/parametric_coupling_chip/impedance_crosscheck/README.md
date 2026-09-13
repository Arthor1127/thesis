# CPW lumped-port impedance cross-check (native PALACE, no SQDMetal)

## The question

SQDMetal's bare-CPW-line driven test (`../driven_sweep/check_line_impedance.py`,
`line_only` mode) backs out `Zin ≈ 33–40 Ω` for a line whose analytic
characteristic impedance is ~50 Ω, and this plateaus under mesh / vacuum-height
/ ground-width refinement. `QC_STATUS.md` still lists this as open (the
"missing ×2 port factor" hypothesis was tested and refuted). Is it a SQDMetal
port-setup artifact, or does PALACE itself get CPW lumped-port impedance wrong?

## What the PALACE docs already say

`palace-src/docs/src/examples/cpw.md` (the developers' own CPW benchmark):

> *"the lumped port excitation exhibits much higher reflection than for wave
> ports. This is expected when using a lumped port to approximate the
> termination of a CPW, and refining the mesh or increasing the order of the
> solution approximation leads to less reflection."* — shown with `p = 2` vs
> `p = 4`.

> *"A multielement uniform lumped port can be used to terminate the CPW by
> connecting the center conductor to the ground plane on each side with
> impedance `Z = 2·Z₀`."*

Their `cpw_lumped_uniform.json` sets `"R": 56.02` (= `Z₀`, the **total** port
impedance, 2 elements) and `"Order": 2`. **SQDMetal's `check_line_impedance.py`
runs `solver_order = 1`.** Prime suspect: the ~40 Ω is the known lumped-port /
CPW discretisation error, worst at order 1.

## The two tests (both LOCAL, both deliberately tiny)

Standing rule is no local FEM. This is a one-off, user-authorised exception for
a single short CPW line. Every launcher `source`s `_guard.sh`, which **requires
`IMPEDANCE_LOCAL_OK=1`**, sets a hard `ulimit -v` (9 GiB → clean `bad_alloc`,
never an OOM-kill), pins threads, caps wall time, and logs `free -m` to
`mem_*.log`. Watch that log — sustained swap growth means `Ctrl-C` now.

```bash
export PALACE_BIN=/home/ruiz/repo/spack/opt/spack/linux-zen3/palace-develop-blpbtxe43ne2or3yo7gou6ycmdwdxuzl/bin/palace
export IMPEDANCE_LOCAL_OK=1
```

### Stage 1 — PALACE's own benchmark, pre-built mesh, zero new geometry

```bash
./stage1_run_palace_cpw.sh        # runs cpw_lumped_uniform at Order 2 and Order 4
python stage1_analyze.py          # Zin from S11 vs the documented 56.02 Ω
```

If PALACE reproduces `Z₀ = 56.02 Ω` (|S11| → −∞) at order 2/4 on its own
matched benchmark, PALACE's lumped-port model is sound and the SQDMetal number
is a SQDMetal issue (order 1, port geometry, mesh, or the R value written).

### Stage 2 — the thesis' own Si CPW cross-section, 2-port through line

```bash
JULIA_PROJECT=/home/ruiz/PalaceForCQED ./stage2_run.sh   # Gmsh.jl mesh + solve, Order 2 & 4
python stage2_analyze.py                                  # exact Z0, eps_eff from (S11,S21)
```

`w = 20 µm`, `s = 11.5 µm`, Si `εr = 11.45`, a 2-fin lumped port (`R = 50`,
total) at **each** end of a straight line — the same measurement as
`../driven_sweep/check_line_impedance.py`. `stage2_analyze.py` runs the exact
uniform-line inversion `(S11, S21) → Z0(f), γ(f)` (NRW-style, no fitting;
verified exact on synthetic data even at a 13-point / 2 GHz sweep) and compares
`Z0` to the analytic Visser value (`analytic_cpw.py`, ported from
`palace-src/examples/cpw/cpw.jl`, matches the Julia to 3 d.p.).

Start with `FREQ_STEP_GHZ=2.0` (13 solves). If `Z0` looks sane, rerun with
`FREQ_STEP_GHZ=0.5` for a tighter number. Analytic target: **Z0 ≈ 50.3 Ω**.

## Files

| file | what |
|---|---|
| `_guard.sh` | safety preamble (`ulimit`, interlock, mem logger) — sourced, not run |
| `analytic_cpw.py` | Visser CPW `Z₀(w,s,h,εr)`, from `cpw.jl` |
| `stage1_run_palace_cpw.sh` / `stage1_analyze.py` | native PALACE CPW benchmark |
| `mesh_si_cpw_open_line.jl` | single open-ended Si CPW line mesh (Gmsh.jl) |
| `gen_mesh.jl` | CLI wrapper for the above |
| `si_cpw_open_line.json` | PALACE driven config template |
| `stage2_run.sh` / `stage2_analyze.py` | Si-line solve + `Z₀` fit |

## Status

Written 2026-09-10. **Not yet run** (awaiting the go-ahead). The Gmsh.jl mesh
script (`mesh_si_cpw_open_line.jl`) has not been executed — expect to iterate on
it once Stage 1 (which needs no new geometry) is confirmed working.
