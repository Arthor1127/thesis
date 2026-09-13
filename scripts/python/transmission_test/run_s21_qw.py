"""
Driven S-parameter (scattering) simulation of the CALIBRATED 5 GHz quarter-wave
resonator, capacitively coupled to a through feedline.

Geometry = batch_test/build_design_quarterwave.py's resonator (open-short,
cpw_width = 20 um, cpw_gap = 12.25 um, total_length = 5.92905 mm - the length
that its elliptic-integral CPW model puts at f_lambda/4 = 5.0 GHz) PLUS a
feedline (LP1 - TL - LP2) that the open end couples to. Lumped ports on the two
launchpads; PROM adaptive fast sweep.

Env: PALACE_BIN, NSLOTS, OUTDIR, FSTART, FEND, NPTS, RES_MIN_UM, OTG1_Y_UM.
"""
import os
import sys
import time
import json

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

OUTDIR     = os.environ.get("OUTDIR", "./s21_qw")
PALACE_BIN = os.environ["PALACE_BIN"]
NUM_CPUS   = int(os.environ.get("NSLOTS", "32"))
os.makedirs(OUTDIR, exist_ok=True)

# --- calibrated resonator (from batch_test/build_design_quarterwave.py) ---
TOTAL_LENGTH_MM = float(os.environ.get("TOTAL_LENGTH_MM", "5.92905"))
CPW_W_UM, CPW_G_UM = 20.0, 12.25
OTG1_Y_UM  = float(os.environ.get("OTG1_Y_UM", "-65"))   # open end; feedline sits at y=0
START_STRAIGHT_UM = 300.0
FILLET_UM = 99.9

# --- sweep ---  eigenmode of the 7.48 mm / 20-12.25 line converged at 3.957 GHz,
# scaling by 7.48345/5.92905 puts this length's lambda/4 at ~4.99 GHz.
FSTART_GHZ = float(os.environ.get("FSTART", "4.0"))
FEND_GHZ   = float(os.environ.get("FEND", "6.0"))
N_POINTS   = int(os.environ.get("NPTS", "2000"))
# PROM=1 (default): Palace's adaptive fast sweep (reduced-order model).
# PROM=0: brute-force full FEM solve at every frequency point -- slow, but the
# ground truth to validate the PROM notch against.
USE_PROM = os.environ.get("PROM", "1") != "0"
ADAPTIVE_TOL, ADAPTIVE_MAX_SAMPLES = 1e-3, 40
RES_MIN_UM, TL_MIN_UM, BG_MAX_UM = float(os.environ.get("RES_MIN_UM", "4.0")), 40.0, 200.0
SOLVER_ORDER = int(os.environ.get("SOLVER_ORDER", "1"))

print(f"L={TOTAL_LENGTH_MM}mm  cpw {CPW_W_UM}/{CPW_G_UM}um  otg1_y={OTG1_Y_UM}um  "
      f"sweep {FSTART_GHZ}-{FEND_GHZ} GHz x{N_POINTS}  res_min={RES_MIN_UM}um", flush=True)


def build():
    from qiskit_metal import designs, Dict
    from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
    from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder
    from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
    from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
    from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround

    design = designs.DesignPlanar({}, overwrite_enabled=True)
    design.chips.main.size.size_x = '4.0mm'
    design.chips.main.size.size_y = '3.0mm'
    design.chips.main.size.size_z = '500um'
    design.chips.main.size.center_x = '0mm'
    design.chips.main.size.center_y = '-0.75mm'
    design.variables['cpw_width'] = f'{CPW_W_UM} um'
    design.variables['cpw_gap'] = f'{CPW_G_UM} um'

    LaunchpadWirebond(design, 'LP1', options=dict(
        chip='main', pos_x='-1.5mm', pos_y='0mm', orientation='360',
        lead_length='30um', pad_height='103um', pad_width='103um', pad_gap='60um'))
    LaunchpadWirebond(design, 'LP2', options=dict(
        chip='main', pos_x='1.5mm', pos_y='0mm', orientation='180',
        lead_length='30um', pad_height='103um', pad_width='103um', pad_gap='60um'))
    RoutePathfinder(design, 'TL', options=dict(
        chip='main', trace_width=f'{CPW_W_UM}um', trace_gap=f'{CPW_G_UM}um', fillet='90um',
        hfss_wire_bonds=True, lead=dict(end_straight='0.1mm'),
        pin_inputs=Dict(start_pin=Dict(component='LP1', pin='tie'),
                        end_pin=Dict(component='LP2', pin='tie'))))

    # open end (couples to the feedline), short end -> lambda/4
    OpenToGround(design, 'otg1', options=dict(
        chip='main', pos_x='-0.2mm', pos_y=f'{OTG1_Y_UM}um', orientation=180,
        width=f'{CPW_W_UM}um', gap=f'{CPW_G_UM}um', termination_gap=f'{CPW_G_UM}um'))
    ShortToGround(design, 'sg1', options=dict(
        chip='main', pos_x='0mm', pos_y='-1.35mm', orientation=-90,
        width=f'{CPW_W_UM}um', gap=f'{CPW_G_UM}um'))

    RouteMeander(design, 'resonator1', Dict(
        trace_width=f'{CPW_W_UM}um', trace_gap=f'{CPW_G_UM}um',
        total_length=f'{TOTAL_LENGTH_MM}mm', hfss_wire_bonds=False,
        fillet=f'{FILLET_UM} um', lead=dict(start_straight=f'{START_STRAIGHT_UM}um'),
        pin_inputs=Dict(start_pin=Dict(component='sg1', pin='short'),
                        end_pin=Dict(component='otg1', pin='open'))))
    return design


t0 = time.time()
design = build()
import matplotlib
matplotlib.use("Agg", force=True)
design.rebuild()
print("components:", list(design.components.keys()), f"  ({time.time()-t0:.1f}s)", flush=True)
r = design.components['resonator1']
print("resonator1 status:", r.status, flush=True)

from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation
driven_sim = PALACE_Driven_Simulation(
    name="s21_qw", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "dielectric_material": "silicon", "solns_to_save": 0,
        "solver_order": SOLVER_ORDER, "solver_tol": 1e-6, "solver_maxits": 500,
        "solver_initial_guess": True, "fillet_resolution": 12,
        "palace_dir": PALACE_BIN, "num_cpus": NUM_CPUS,
    })
driven_sim.add_metallic(1)
driven_sim.add_ground_plane()
driven_sim.create_port_CPW_on_Launcher("LP1", 30e-6)
driven_sim.create_port_CPW_on_Launcher("LP2", 30e-6)
driven_sim.set_port_excitation(1)
driven_sim.fine_mesh_components(["TL", "LP1", "LP2"], min_size=TL_MIN_UM*1e-6, max_size=BG_MAX_UM*1e-6,
                                taper_dist_min=2*TL_MIN_UM*1e-6, taper_dist_max=5*TL_MIN_UM*1e-6, metals_only=False)
driven_sim.fine_mesh_components(["resonator1", "otg1", "sg1"], min_size=RES_MIN_UM*1e-6, max_size=BG_MAX_UM*1e-6,
                                taper_dist_min=30e-6, taper_dist_max=100e-6, metals_only=False)
# coupling region: a MODEST refinement of the ground strip between the feedline
# (y=0) and the resonator top row (y~OTG1_Y_UM), only where they run closest
# (|x| < 0.6 mm). Keep it gentle - a 1.5 um box over the whole strip took the
# mesh to 17 M elements and segfaulted SuperLU.
_cbox = float(os.environ.get("CPL_BOX_UM", "2.5"))
if _cbox > 0:
    driven_sim.fine_mesh_in_rectangle(-0.6e-3, OTG1_Y_UM*1e-6 - 20e-6, 0.6e-3, -15e-6,
                                      min_size=_cbox*1e-6, max_size=BG_MAX_UM*1e-6,
                                      taper_dist_min=20e-6, taper_dist_max=70e-6)

fs, fe = FSTART_GHZ*1e9, FEND_GHZ*1e9
driven_sim.set_freq_values(freq_start=fs, freq_end=fe, freq_step=(fe-fs)/max(N_POINTS-1, 1))

tp = time.time()
driven_sim.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

cfg_path = driven_sim._sim_config
with open(cfg_path) as f:
    cfg = json.load(f)
if USE_PROM:
    cfg["Solver"]["Driven"]["AdaptiveTol"] = ADAPTIVE_TOL
    cfg["Solver"]["Driven"]["AdaptiveMaxSamples"] = ADAPTIVE_MAX_SAMPLES
else:
    cfg["Solver"]["Driven"].pop("AdaptiveTol", None)
    cfg["Solver"]["Driven"].pop("AdaptiveMaxSamples", None)
with open(cfg_path, "w") as f:
    json.dump(cfg, f, indent=2)
print(f"PROM={'on' if USE_PROM else 'OFF (full solve every point)'}:",
      cfg["Solver"]["Driven"], flush=True)

ts = time.time()
driven_sim.run()
print(f"SOLVE: {time.time()-ts:.1f} s", flush=True)

s = driven_sim.retrieve_data()
if s is None:
    raise SystemExit("retrieve_data() returned None - the Palace solve failed "
                     "(check run.log for a segfault / SuperLU MemType error; "
                     "likely the mesh is too large for the direct solver).")
freqs = np.asarray(s["freqs"]); S21 = np.asarray(s["S21"]); S11 = np.asarray(s.get("S11"))
_tag = "prom" if USE_PROM else "full"
np.savez(f"{OUTDIR}/s21_qw_{_tag}.npz", freq_hz=freqs, S21=S21, S11=S11,
         total_length_mm=TOTAL_LENGTH_MM, otg1_y_um=OTG1_Y_UM, prom=USE_PROM)
s21_db = 20*np.log10(np.abs(S21))
i = int(np.argmin(s21_db))
print(f"|S21| min {s21_db[i]:.3f} dB @ {freqs[i]/1e9:.5f} GHz", flush=True)
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
