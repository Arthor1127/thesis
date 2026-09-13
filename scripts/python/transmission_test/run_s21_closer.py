"""
Driven (S21) run with the resonator moved CLOSER to the feedline, to deepen the
(currently 0.23 dB) coupling notch.

Self-contained: builds its own geometry (a parametrized copy of build_design.py's
body) so it doesn't depend on the WIP build_design.py signature. Uses the PROM
fast sweep (AdaptiveTol) and order 1.

Env: PALACE_BIN, NSLOTS, OUTDIR, plus the geometry knobs OTG1_Y_UM, STG2_Y_MM,
START_STRAIGHT_UM, and the sweep FSTART/FEND/NPTS.
"""
import os
import sys
import time
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

OUTDIR     = os.environ.get("OUTDIR", "./s21_closer")
PALACE_BIN = os.environ["PALACE_BIN"]
NUM_CPUS   = int(os.environ.get("NSLOTS", "32"))
os.makedirs(OUTDIR, exist_ok=True)

# --- geometry knobs (defaults = current build_design.py; override to move closer) ---
OTG1_Y_UM        = float(os.environ.get("OTG1_Y_UM", "-40"))     # open end: -40 -> up toward the feedline (y=0)
STG2_Y_MM        = float(os.environ.get("STG2_Y_MM", "-1.25"))   # short end: -1.25 -> up
START_STRAIGHT_UM = float(os.environ.get("START_STRAIGHT_UM", "10"))  # open-end lead straight (longer = more coupling length)
TOTAL_LENGTH_MM  = 7.48345

FSTART_GHZ = float(os.environ.get("FSTART", "3.2"))
FEND_GHZ   = float(os.environ.get("FEND", "4.4"))
N_POINTS   = int(os.environ.get("NPTS", "3001"))
ADAPTIVE_TOL, ADAPTIVE_MAX_SAMPLES = 1e-3, 30
RES_MIN_UM, TL_MIN_UM, BG_MAX_UM = 6.0, 40.0, 200.0
SOLVER_ORDER = 1

print(f"otg1_y={OTG1_Y_UM}um  stg2_y={STG2_Y_MM}mm  start_straight={START_STRAIGHT_UM}um  "
      f"sweep {FSTART_GHZ}-{FEND_GHZ} GHz x{N_POINTS}", flush=True)


def build_closer():
    from qiskit_metal import designs, Dict
    from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
    from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder
    from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
    from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
    from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround

    design = designs.DesignPlanar({}, overwrite_enabled=True)
    design.chips.main.size.size_x = '4.0mm'
    design.chips.main.size.size_y = '3mm'
    design.chips.main.size.size_z = '500um'
    design.chips.main.size.center_x = '0mm'
    design.chips.main.size.center_y = '-0.5mm'
    design.variables['cpw_width'] = '10 um'
    design.variables['cpw_gap'] = '6 um'

    LaunchpadWirebond(design, 'LP1', options=dict(
        chip='main', pos_x='-1mm', pos_y='0mm', orientation='360',
        lead_length='30um', pad_height='103um', pad_width='103um', pad_gap='60um'))
    LaunchpadWirebond(design, 'LP2', options=dict(
        chip='main', pos_x='1mm', pos_y='0mm', orientation='180',
        lead_length='30um', pad_height='103um', pad_width='103um', pad_gap='60um'))
    RoutePathfinder(design, 'TL', options=dict(
        chip='main', trace_width='10um', trace_gap='6um', fillet='90um',
        hfss_wire_bonds=True, lead=dict(end_straight='0.1mm'),
        pin_inputs=Dict(start_pin=Dict(component='LP1', pin='tie'),
                        end_pin=Dict(component='LP2', pin='tie'))))

    OpenToGround(design, 'otg1', options=dict(
        chip='main', pos_x='-0.2mm', pos_y=f'{OTG1_Y_UM}um', orientation=180))
    ShortToGround(design, 'stg2', options=dict(
        chip='main', pos_x='0mm', pos_y=f'{STG2_Y_MM}mm', orientation=0))

    RouteMeander(design, 'resonator1', Dict(
        trace_width='10um', trace_gap='6um', total_length=f'{TOTAL_LENGTH_MM}mm',
        hfss_wire_bonds=False, fillet='99.9 um',
        lead=dict(start_straight=f'{START_STRAIGHT_UM}um'),
        meander=Dict(spacing="200um", assymetry="1mm"),
        pin_inputs=Dict(start_pin=Dict(component='otg1', pin='open'),
                        end_pin=Dict(component='stg2', pin='short'))))
    return design


t0 = time.time()
design = build_closer()
import matplotlib
matplotlib.use("Agg", force=True)
design.rebuild()
print("components:", list(design.components.keys()), f"  ({time.time()-t0:.1f}s)", flush=True)

# report open-end proximity to the feedline (feedline gap edge ~ y = -11 um for cpw_gap=6, width=10)
print(f"open end at y={OTG1_Y_UM}um; feedline outer gap edge ~ y=-11um -> ~{abs(OTG1_Y_UM)-11:.0f}um of ground between", flush=True)

from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation
driven_sim = PALACE_Driven_Simulation(
    name="s21_closer", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "dielectric_material": "silicon", "solns_to_save": 0,
        "solver_order": SOLVER_ORDER, "solver_tol": 1e-6, "solver_maxits": 500,
        "solver_initial_guess": True, "fillet_resolution": 12,
        "palace_dir": PALACE_BIN, "num_cpus": NUM_CPUS,
    })
driven_sim.add_metallic(1)
driven_sim.add_ground_plane()
driven_sim.create_port_CPW_on_Launcher("LP1", 20e-6)
driven_sim.create_port_CPW_on_Launcher("LP2", 20e-6)
driven_sim.set_port_excitation(1)
driven_sim.fine_mesh_components(["resonator1"], min_size=RES_MIN_UM*1e-6, max_size=BG_MAX_UM*1e-6,
                                taper_dist_min=3*RES_MIN_UM*1e-6, taper_dist_max=8*RES_MIN_UM*1e-6, metals_only=False)
driven_sim.fine_mesh_components(["TL", "LP1", "LP2"], min_size=TL_MIN_UM*1e-6, max_size=BG_MAX_UM*1e-6,
                                taper_dist_min=2*TL_MIN_UM*1e-6, taper_dist_max=5*TL_MIN_UM*1e-6, metals_only=False)

fs, fe = FSTART_GHZ*1e9, FEND_GHZ*1e9
driven_sim.set_freq_values(freq_start=fs, freq_end=fe, freq_step=(fe-fs)/max(N_POINTS-1, 1))

tp = time.time()
driven_sim.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

cfg_path = driven_sim._sim_config
with open(cfg_path) as f:
    cfg = json.load(f)
cfg["Solver"]["Driven"]["AdaptiveTol"] = ADAPTIVE_TOL
cfg["Solver"]["Driven"]["AdaptiveMaxSamples"] = ADAPTIVE_MAX_SAMPLES
with open(cfg_path, "w") as f:
    json.dump(cfg, f, indent=2)
print("PROM:", cfg["Solver"]["Driven"], flush=True)

ts = time.time()
driven_sim.run()
print(f"SOLVE: {time.time()-ts:.1f} s", flush=True)

s = driven_sim.retrieve_data()
freqs = np.asarray(s["freqs"]); S21 = np.asarray(s["S21"]); S11 = np.asarray(s.get("S11"))
np.savez(f"{OUTDIR}/s21_closer.npz", freq_hz=freqs, S21=S21, S11=S11,
         otg1_y_um=OTG1_Y_UM, stg2_y_mm=STG2_Y_MM, start_straight_um=START_STRAIGHT_UM)
s21_db = 20*np.log10(np.abs(S21))
i = int(np.argmin(s21_db))
print(f"|S21| min {s21_db[i]:.3f} dB @ {freqs[i]/1e9:.5f} GHz   (baseline was -0.23 dB @ 3.8265 GHz)", flush=True)
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
