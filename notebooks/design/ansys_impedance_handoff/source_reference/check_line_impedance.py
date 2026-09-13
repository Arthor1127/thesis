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

CPW_WIDTH_UM = 20.0
CPW_GAP_UM = 11.5
CHIP_HALF_SPAN_UM = 1000.0

design = designs.DesignPlanar({}, overwrite_enabled=True)
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

driven = PALACE_Driven_Simulation(
    name=f"line_{MODE}", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "dielectric_material": "silicon", "solns_to_save": 0,
        "solver_order": 1, "solver_tol": 1e-6, "solver_maxits": 500,
        "solver_initial_guess": True, "fillet_resolution": 12,
        "palace_dir": PALACE_BIN, "num_cpus": NCPU,
    })
driven.add_metallic(1)
driven.add_ground_plane()

if MODE == "with_launchpads":
    driven.create_port_CPW_on_Launcher("LP1", 30e-6)
    driven.create_port_CPW_on_Launcher("LP2", 30e-6)
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
    driven.create_port_CPW_on_Route("TL", pin_name="start", len_launch=30e-6)
    driven.create_port_CPW_on_Route("TL", pin_name="end", len_launch=30e-6)
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
MIN_SIZE_FACTOR = float(os.environ.get("MIN_SIZE_FACTOR", "0.5"))
Z_PROP_TOP = float(os.environ.get("Z_PROP_TOP", "1.0"))
X_PROP = float(os.environ.get("X_PROP", "0.1"))
Y_PROP = float(os.environ.get("Y_PROP", "0.1"))
lateral_min_um = MIN_SIZE_FACTOR * min(CPW_WIDTH_UM, CPW_GAP_UM)
driven.fine_mesh_components(fine_comp_names, min_size=lateral_min_um * 1e-6,
                            max_size=100e-6, taper_dist_min=16e-6, taper_dist_max=40e-6,
                            metals_only=False)
driven.set_zBoundary_as_proportion(Z_PROP_TOP, 0.0)
# X_PROP/Y_PROP: the mesh x2/z x2 convergence sweep (see
# submit_line_convergence.sh) plateaued at |Zin|=40.483 Ohm regardless of
# further mesh refinement OR vacuum height -- neither is the bottleneck.
# Lateral ground-plane width is the remaining untested truncation: defaults
# (x_prop=y_prop=0.1) give only 200um of ground beyond the chip edge in y
# (total half-width ~1.2mm from the trace) and 480um in x.
driven.set_xBoundary_as_proportion(X_PROP)
driven.set_yBoundary_as_proportion(Y_PROP)
print(f"mesh: lateral_min={lateral_min_um:.3f}um (factor={MIN_SIZE_FACTOR})  "
     f"z_prop_top={Z_PROP_TOP} (vacuum height above chip = {Z_PROP_TOP*500}um)  "
     f"x_prop={X_PROP} y_prop={Y_PROP}", flush=True)

# No resonator -> no resonance to find -> no need for PROM/dense adaptive
# sampling. A short, coarse frequency list across the same band used for the
# resonator sweep is enough to confirm Zin is flat and near 50 Ohm everywhere.
FREQS_GHZ = [6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
fs, fe = FREQS_GHZ[0] * 1e9, FREQS_GHZ[-1] * 1e9
driven.set_freq_values(freq_start=fs, freq_end=fe, freq_step=(fe - fs) / (len(FREQS_GHZ) - 1))

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

Z0_ref = 50.0
Zin = Z0_ref * (1 + S11) / (1 - S11)
row = {
    "mode": MODE, "tag": TAG, "min_size_factor": MIN_SIZE_FACTOR,
    "lateral_min_um": lateral_min_um, "z_prop_top": Z_PROP_TOP,
    "x_prop": X_PROP, "y_prop": Y_PROP,
    "freqs_ghz": (freqs / 1e9).tolist(),
    "S11_db": (20 * np.log10(np.abs(S11))).tolist(),
    "Zin_real": Zin.real.tolist(), "Zin_imag": Zin.imag.tolist(),
    "mean_Zin_abs": float(np.mean(np.abs(Zin))),
    "mean_Zin_real": float(np.mean(Zin.real)), "mean_Zin_imag": float(np.mean(Zin.imag)),
    "peak_mem_per_rank_total_gb": peak_mem_per_rank_gb, "peak_mem_per_node_gb": peak_mem_total_gb,
    "solve_s": solve_s, "wall_s": time.time() - t0,
}
with open(f"{OUTDIR}/results/{TAG}.json", "w") as f:
    json.dump(row, f, indent=2)
print(json.dumps({k: v for k, v in row.items() if not isinstance(v, list)}, indent=2), flush=True)
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
