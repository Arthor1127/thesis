"""Bare-CPW-line impedance check, with NO resonator/coupler in the design at
all: isolates whether Zin!=50Ohm (found for the full cs14_cl300 design: Zin=
20.9+5.5j Ohm away from resonance) comes from the CoupledLineTee coupler, the
LaunchpadWirebond tapers, or the plain CPW cross-section itself.

Two modes, run as SEPARATE jobs so each isolates one thing:

  --mode with_launchpads : LP1 -- straight CPW route -- LP2 (real
      LaunchpadWirebond tapers, ports via create_port_CPW_on_Launcher, same
      convention as run_qc_sim.py/readout_resonators.py). Isolates:
      launchpad-taper + bare-line contribution, with NO coupler.

  --mode line_only : minimal OpenToGround stubs at each end (just to satisfy
      RoutePathfinder's pin requirement) -- straight CPW route -- ports placed
      via create_port_CPW_via_edge_point a few hundred um IN from each stub,
      on the plain trace itself. No LaunchpadWirebond, no coupler. Isolates:
      JUST the CPW cross-section's own Z0.

No resonator means no PROM adaptive resonance-finding is needed -- a matched
50Ohm line should show S11~0 and S21~1 at every frequency, so this uses a
short, coarse, cheap frequency list rather than run_qc_sim.py's dense window.

Env: MODE (with_launchpads|line_only), PALACE_BIN, [NSLOTS=16], [OUTDIR=./out_line]
Writes: {OUTDIR}/results/{MODE}.json, {OUTDIR}/results/{MODE}_s21.npz
Cluster job only -- run via submit_line_impedance.sh.
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

MODE = os.environ["MODE"]
assert MODE in ("with_launchpads", "line_only")
OUTDIR = os.environ.get("OUTDIR", "./out_line")
PALACE_BIN = os.environ["PALACE_BIN"]
NCPU = int(os.environ.get("NSLOTS", "16"))
TAG = os.environ.get("TAG", MODE)
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

print(f"[{MODE}] ncpu={NCPU}", flush=True)

t0 = time.time()
from qiskit_metal import designs, Dict  # noqa: F401
import matplotlib
matplotlib.use("Agg", force=True)

CPW_WIDTH_UM = float(os.environ.get("CPW_WIDTH_UM", "20.0"))
CPW_GAP_UM = float(os.environ.get("CPW_GAP_UM", "11.5"))
CHIP_HALF_SPAN_UM = 1000.0

# MATCH_JULIA: when set, line_only's chip geometry is derived EXACTLY like
# impedance_crosscheck/mesh_si_cpw_line.jl -- substrate/chip un-padded in x/y
# (edges flush with the ports and the ground-trace-ground pattern, size_y =
# 2*GROUND_UM+2*CPW_GAP_UM+CPW_WIDTH_UM, same as Julia's y5), same substrate
# thickness, ports at the chip edges (CHIP_HALF_SPAN_UM = LENGTH_UM/2), and
# the same 400um air pad in x/y/z (including below the substrate -- Julia's
# substrate floats in open air on both z faces, no grounded backplane).
# Previously this script used an arbitrary 4.8mm x 2mm x 500um chip with
# ports 2000um apart sitting deep inside it -- not the same structure at all.
MATCH_JULIA = os.environ.get("MATCH_JULIA", "0") == "1"
GROUND_UM = float(os.environ.get("GROUND_UM", "250.0"))
LENGTH_UM = float(os.environ.get("LENGTH_UM", "4000.0"))
SUBSTRATE_H_UM = float(os.environ.get("SUBSTRATE_H_UM", "525.0"))
AIR_PAD_UM = float(os.environ.get("AIR_PAD_UM", "400.0"))
Y5_UM = 2 * GROUND_UM + 2 * CPW_GAP_UM + CPW_WIDTH_UM

design = designs.DesignPlanar({}, overwrite_enabled=True)
if MATCH_JULIA:
    CHIP_HALF_SPAN_UM = LENGTH_UM / 2.0
    design.chips.main.size.size_x = f'{LENGTH_UM}um'
    design.chips.main.size.size_y = f'{Y5_UM}um'
    design.chips.main.size.size_z = f'{SUBSTRATE_H_UM}um'
else:
    design.chips.main.size.size_x = '4.8mm'
    design.chips.main.size.size_y = '2mm'
    design.chips.main.size.size_z = '500um'
design.chips.main.size.center_x = '0mm'
design.chips.main.size.center_y = '0mm'
design.variables['cpw_width'] = f'{CPW_WIDTH_UM}um'
design.variables['cpw_gap'] = f'{CPW_GAP_UM}um'

if MODE == "with_launchpads":
    from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
    from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder

    x1, x2 = f'-{CHIP_HALF_SPAN_UM}um', f'{CHIP_HALF_SPAN_UM}um'
    LaunchpadWirebond(design, 'LP1', options=dict(
        chip='main', pos_x=x1, pos_y='0um', orientation='360',
        lead_length='30um', pad_height='400um', pad_width='250um',
        pad_gap='150um', taper_height='350um'))
    LaunchpadWirebond(design, 'LP2', options=dict(
        chip='main', pos_x=x2, pos_y='0um', orientation='180',
        lead_length='30um', pad_height='400um', pad_width='250um',
        pad_gap='150um', taper_height='350um'))
    RoutePathfinder(design, 'TL', options=dict(
        chip='main', trace_width='cpw_width', trace_gap='cpw_gap', fillet='50um',
        hfss_wire_bonds=True,
        pin_inputs=Dict(start_pin=Dict(component='LP1', pin='tie'),
                        end_pin=Dict(component='LP2', pin='tie'))))
    design.rebuild()
    print(f"geometry built ({time.time()-t0:.1f}s): LP1--TL--LP2, no coupler", flush=True)
else:
    from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
    from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder

    x1, x2 = f'-{CHIP_HALF_SPAN_UM}um', f'{CHIP_HALF_SPAN_UM}um'
    OpenToGround(design, 'OTG1', options=dict(
        chip='main', pos_x=x1, pos_y='0um', orientation='180', width='cpw_width'))
    OpenToGround(design, 'OTG2', options=dict(
        chip='main', pos_x=x2, pos_y='0um', orientation='0', width='cpw_width'))
    RoutePathfinder(design, 'TL', options=dict(
        chip='main', trace_width='cpw_width', trace_gap='cpw_gap', fillet='50um',
        pin_inputs=Dict(start_pin=Dict(component='OTG1', pin='open'),
                        end_pin=Dict(component='OTG2', pin='open'))))
    design.rebuild()
    print(f"geometry built ({time.time()-t0:.1f}s): OTG1--TL--OTG2, no launchpads, "
          f"no coupler", flush=True)

from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation

# MIN_SIZE_FACTOR computed here (before constructing `driven`) so
# gmsh_dist_func_discretisation can be derived from it below.
MIN_SIZE_FACTOR = float(os.environ.get("MIN_SIZE_FACTOR",
    str(0.75 * CPW_WIDTH_UM / min(CPW_WIDTH_UM, CPW_GAP_UM)) if MATCH_JULIA else "0.5"))
lateral_min_um = MIN_SIZE_FACTOR * min(CPW_WIDTH_UM, CPW_GAP_UM)

# GMSH_DIST_SAMPLING: Julia's Distance field sets Sampling=ceil(L/l_trace),
# scaling the reference-geometry discretization to the target near-line
# element size (~500+ samples for this 4mm line). SQDMetal's
# GMSH_Mesh_Builder.py leaves this at a flat default of 150 regardless of
# geometry size -- for a 4mm-long "TL" component that's ~27um between
# samples, far coarser than the few-um target element size, which plausibly
# explains the blockier (vs Julia's smooth cone) fine-mesh region observed.
# Matching Julia's scaling convention under MATCH_JULIA.
GMSH_DIST_SAMPLING = int(float(os.environ.get("GMSH_DIST_SAMPLING",
    str(LENGTH_UM / lateral_min_um) if MATCH_JULIA else "150")))

driven = PALACE_Driven_Simulation(
    name=f"line_{MODE}", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "dielectric_material": "silicon", "solns_to_save": 0,
        "solver_order": int(os.environ.get("SOLVER_ORDER", "1")), "solver_tol": 1e-6, "solver_maxits": 500,
        "solver_initial_guess": True, "fillet_resolution": 12,
        "palace_dir": PALACE_BIN, "num_cpus": NCPU,
        "gmsh_dist_func_discretisation": GMSH_DIST_SAMPLING,
    })
driven.add_metallic(1)
driven.add_ground_plane()

# FARFIELD_TYPE: set_farfield() is never called anywhere in this pipeline by
# default, and Frequency_Driven_Simulation.__init__ sets self._ff_type = {}
# -- confirmed via the actual generated config.json (grep -c '"Absorbing"'
# returns 0) that every outer air-box face (all 6, not just the z=0 floor)
# defaults to PEC, only 200-500um from the CPW (x_prop/y_prop=0.1,
# z_prop_top=1.0). That's a fully shielded/boxed CPW, not the open line the
# Visser analytic Z0 assumes -- a much bigger deviation than port geometry or
# ground-plane width. The Julia raw-Gmsh crosscheck uses an explicit
# Absorbing boundary; matching that here.
FARFIELD_TYPE = os.environ.get("FARFIELD_TYPE", "absorbing")
driven.set_farfield(FARFIELD_TYPE)

# PORT_IMPEDANCE_R: root-cause test. PALACE's own CPW example docs
# (docs/src/examples/cpw.md) state a multielement uniform lumped port (the
# "fins to ground from the central CPW stripline" pattern
# create_port_CPW_on_Launcher/on_Route both use -- see their docstrings)
# must have EACH element set to Z=2*Z0 to realize a net port impedance of
# Z0, because the two fin elements are electrically in parallel (2Z0 || 2Z0
# = Z0). Checked SQDMetal's own source (Model.py:1708-1709): it writes
# leDict['R'] = impedance_R straight through with NO x2 compensation for
# this. Default here (50) reproduces every earlier run exactly; passing 100
# tests whether that missing factor of 2 is the actual cause of the
# Zin~40 Ohm (vs analytic ~50 Ohm) discrepancy chased all session.
PORT_IMPEDANCE_R = float(os.environ.get("PORT_IMPEDANCE_R", "50"))

if MODE == "with_launchpads":
    driven.create_port_CPW_on_Launcher("LP1", 30e-6, impedance_R=PORT_IMPEDANCE_R)
    driven.create_port_CPW_on_Launcher("LP2", 30e-6, impedance_R=PORT_IMPEDANCE_R)
    fine_comp_names = ["TL", "LP1", "LP2"]
else:
    # create_port_CPW_via_edge_point needs an actual open-circuit END of a
    # conductor to detect a clean 2-point perpendicular edge -- it can't be
    # used mid-trace (tried inset 200um from the OTG stubs first; it threw
    # "Seems that a CPW couldn't be formed" because there's no discontinuity
    # there at all, just a continuous line). create_port_CPW_on_Route attaches
    # directly to the Route's own start/end pin instead, which IS a real end
    # -- physically the same location as the OTG stub, but OTG's own geometry
    # is just a thin gap-etch cap (no taper), so this still isolates the bare
    # line with no launchpad-taper contribution.
    # LEN_LAUNCH_UM: create_port_CPW_on_Route's own docstring says this should
    # be kept THIN w.r.t. the CPW gap -- the old default (30um) is 5x the
    # 6um gap used by CPW_GAP_UM=6 tests, the opposite of that guidance, and
    # plausibly the source of the large spurious Im(Zin) this mode shows that
    # the geometry-matched Julia crosscheck (fin length = gap width, per
    # Palace's own cpw.jl convention) does not.
    # Default fin length = gap width, matching Palace's own documented square
    # fin convention (and Julia's port fins) when MATCH_JULIA; else keep the
    # old 30um default for other exploratory uses.
    LEN_LAUNCH_UM = float(os.environ.get("LEN_LAUNCH_UM", str(CPW_GAP_UM) if MATCH_JULIA else "30.0"))
    driven.create_port_CPW_on_Route("TL", pin_name="start", len_launch=LEN_LAUNCH_UM * 1e-6,
                                    impedance_R=PORT_IMPEDANCE_R)
    driven.create_port_CPW_on_Route("TL", pin_name="end", len_launch=LEN_LAUNCH_UM * 1e-6,
                                    impedance_R=PORT_IMPEDANCE_R)
    fine_comp_names = ["TL", "OTG1", "OTG2"]
driven.set_port_excitation(1)

# MIN_SIZE_FACTOR/Z_PROP_TOP let this script be used as a convergence-sweep
# probe (see submit_line_impedance_convergence.sh): default reproduces the
# original with_launchpads/line_only runs exactly (factor=0.5 -> 5.75um,
# z_prop_top=1.0 -> SQDMetal's own PALACE_Model_Base default, 500um of vacuum
# above the chip). Investigating why the bare line's simulated Zin (~33 Ohm)
# disagrees with both DesignHelper's and SQDMetal's own CpwParams analytic
# formula (both ~50 Ohm, and CpwParams is explicitly a grounded-backplane-
# aware formula, so the z_prop_bottom=0 PEC floor at 500um below the trace
# isn't the culprit -- that's already accounted for and predicted negligible).
# (MIN_SIZE_FACTOR/lateral_min_um computed earlier, before constructing `driven`.)
if MATCH_JULIA:
    # Same 400um pad in every direction as Julia, expressed as the proportion
    # SQDMetal's API wants (of chip size_x/size_y/size_z respectively) --
    # AND a nonzero Z_PROP_BOTTOM, since Julia's substrate floats in open air
    # on both z faces (no grounded backplane), unlike the z_prop_bottom=0.0
    # default used everywhere else in this script.
    Z_PROP_TOP = float(os.environ.get("Z_PROP_TOP", str(AIR_PAD_UM / SUBSTRATE_H_UM)))
    Z_PROP_BOTTOM = float(os.environ.get("Z_PROP_BOTTOM", str(AIR_PAD_UM / SUBSTRATE_H_UM)))
    X_PROP = float(os.environ.get("X_PROP", str(AIR_PAD_UM / LENGTH_UM)))
    Y_PROP = float(os.environ.get("Y_PROP", str(AIR_PAD_UM / Y5_UM)))
else:
    Z_PROP_TOP = float(os.environ.get("Z_PROP_TOP", "1.0"))
    Z_PROP_BOTTOM = float(os.environ.get("Z_PROP_BOTTOM", "0.0"))
    X_PROP = float(os.environ.get("X_PROP", "0.1"))
    Y_PROP = float(os.environ.get("Y_PROP", "0.1"))
# MAX_SIZE_UM/TAPER_DIST_*_UM: Julia's own sizing field (mesh_si_cpw_line.jl)
# tapers from l_trace to l_farfield=1.0*H*2^-refn over DistMin=W to
# DistMax=0.9*(H+2*pad) -- i.e. a >1mm-wide transition, not the 24um cliff
# (16um to 40um) hardcoded here before. A sharp size jump like that produces
# poorly-shaped transition elements and can bias a near-field-sensitive
# quantity like Z0. Matching Julia's gradual taper under MATCH_JULIA.
if MATCH_JULIA:
    MAX_SIZE_UM = float(os.environ.get("MAX_SIZE_UM", str(1.0 * SUBSTRATE_H_UM)))
    TAPER_DIST_MIN_UM = float(os.environ.get("TAPER_DIST_MIN_UM", str(CPW_WIDTH_UM)))
    TAPER_DIST_MAX_UM = float(os.environ.get("TAPER_DIST_MAX_UM", str(0.9 * (SUBSTRATE_H_UM + 2 * AIR_PAD_UM))))
else:
    MAX_SIZE_UM = float(os.environ.get("MAX_SIZE_UM", "100.0"))
    TAPER_DIST_MIN_UM = float(os.environ.get("TAPER_DIST_MIN_UM", "16.0"))
    TAPER_DIST_MAX_UM = float(os.environ.get("TAPER_DIST_MAX_UM", "40.0"))
driven.fine_mesh_components(fine_comp_names, min_size=lateral_min_um * 1e-6,
                            max_size=MAX_SIZE_UM * 1e-6, taper_dist_min=TAPER_DIST_MIN_UM * 1e-6,
                            taper_dist_max=TAPER_DIST_MAX_UM * 1e-6, metals_only=False)
driven.set_zBoundary_as_proportion(Z_PROP_TOP, Z_PROP_BOTTOM)
driven.set_xBoundary_as_proportion(X_PROP)
driven.set_yBoundary_as_proportion(Y_PROP)
print(f"mesh: lateral_min={lateral_min_um:.3f}um (factor={MIN_SIZE_FACTOR})  "
     f"z_prop_top={Z_PROP_TOP} (vacuum height above chip = {Z_PROP_TOP*500}um)  "
     f"x_prop={X_PROP} y_prop={Y_PROP}", flush=True)

# No resonator -> no resonance to find -> no need for PROM/dense adaptive
# sampling. A short, coarse frequency list across the same band used for the
# resonator sweep is enough to confirm Zin is flat and near 50 Ohm everywhere.
FREQS_GHZ = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0] \
    if MATCH_JULIA else [6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
fs, fe = FREQS_GHZ[0] * 1e9, FREQS_GHZ[-1] * 1e9
driven.set_freq_values(freq_start=fs, freq_end=fe, freq_step=(fe - fs) / (len(FREQS_GHZ) - 1))

# Solution-based AMR, same env-gated pattern as run_qc_sim.py.
AMR_ITERS = int(os.environ.get("AMR_ITERS", "0"))
if AMR_ITERS > 0:
    AMR_TOL = float(os.environ.get("AMR_TOL", "1e-2"))
    AMR_MAX_DOFS = int(os.environ.get("AMR_MAX_DOFS", "1500000"))
    AMR_DORFLER = float(os.environ.get("AMR_DORFLER", "0.7"))
    driven.enable_mesh_refinement(num_iterations=AMR_ITERS, max_DoFs=AMR_MAX_DOFS,
                                  tolerance=AMR_TOL, Dorfler_marking_fraction=AMR_DORFLER,
                                  save_iterations_data=True, save_iterations_mesh=False)
    print(f"Solution-based AMR ENABLED: iters={AMR_ITERS} tol={AMR_TOL:.1e} "
         f"max_DoFs={AMR_MAX_DOFS} Dorfler={AMR_DORFLER}", flush=True)

tp = time.time()
driven.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

cfg_path = driven._sim_config
with open(cfg_path) as f:
    cfg = json.load(f)
cfg["Solver"]["Driven"]["AdaptiveTol"] = 1e-3
cfg["Solver"]["Driven"]["AdaptiveMaxSamples"] = 10  # no resonance to resolve
with open(cfg_path, "w") as f:
    json.dump(cfg, f, indent=2)

ts = time.time()
driven.run()
solve_s = time.time() - ts
print(f"SOLVE: {solve_s:.1f} s", flush=True)

s = driven.retrieve_data()
if s is None:
    raise SystemExit("retrieve_data() returned None -- the Palace solve failed.")

import re
peak_mem_per_rank_gb = peak_mem_total_gb = None
try:
    log_path = os.path.join(driven._output_data_dir, "out.log")
    with open(log_path) as f:
        log_text = f.read()
    m = re.search(r"Estimated peak per-rank memory usage is:.*?Total ([\d.]+)G", log_text)
    if m: peak_mem_per_rank_gb = float(m.group(1))
    m = re.search(r"Estimated peak per-node memory usage is:.*?Total ([\d.]+)G", log_text)
    if m: peak_mem_total_gb = float(m.group(1))
except Exception as e:
    print(f"could not parse peak memory: {e}", flush=True)

import numpy as np
freqs = np.asarray(s["freqs"]); S11 = np.asarray(s["S11"]); S21 = np.asarray(s["S21"])
np.savez(f"{OUTDIR}/results/{TAG}_s21.npz", freq_hz=freqs, S11=S11, S21=S21)

# S11 as returned by retrieve_data() is in PALACE's OWN native reference
# impedance for however the port was configured (i.e. PORT_IMPEDANCE_R, not
# necessarily 50) -- S11_db below is reference-impedance-agnostic (just
# |S11| in dB) and is the most direct match-quality check regardless of the
# x2-factor question. Zin_* below assumes the reference IS PORT_IMPEDANCE_R
# (the value actually passed to create_port_CPW_on_Route/on_Launcher), not a
# hardcoded 50 -- important now that this script can set PORT_IMPEDANCE_R=100.
Z0_ref = PORT_IMPEDANCE_R
Zin = Z0_ref * (1 + S11) / (1 - S11)
row = {
    "mode": MODE, "tag": TAG, "min_size_factor": MIN_SIZE_FACTOR,
    "lateral_min_um": lateral_min_um, "z_prop_top": Z_PROP_TOP,
    "x_prop": X_PROP, "y_prop": Y_PROP, "port_impedance_R": PORT_IMPEDANCE_R,
    "freqs_ghz": (freqs / 1e9).tolist(),
    "S11_db": (20 * np.log10(np.abs(S11))).tolist(),
    "Zin_real": Zin.real.tolist(), "Zin_imag": Zin.imag.tolist(),
    "mean_Zin_abs": float(np.mean(np.abs(Zin))),
    "mean_Zin_real": float(np.mean(Zin.real)), "mean_Zin_imag": float(np.mean(Zin.imag)),
    "mean_S11_db": float(np.mean(20 * np.log10(np.abs(S11)))),
    "peak_mem_per_rank_total_gb": peak_mem_per_rank_gb, "peak_mem_per_node_gb": peak_mem_total_gb,
    "solve_s": solve_s, "wall_s": time.time() - t0,
}
with open(f"{OUTDIR}/results/{TAG}.json", "w") as f:
    json.dump(row, f, indent=2)
print(json.dumps({k: v for k, v in row.items() if not isinstance(v, list)}, indent=2), flush=True)
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
