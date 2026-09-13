"""
Eigenmode run of the BARE resonator (no feedline / launchpads / transmission
line) -> the intrinsic mode frequency + internal Q (set by the substrate
LossTan and the metal/interface participations).

Self-contained geometry so it doesn't depend on build_design.py's signature.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd

OUTDIR     = os.environ.get("OUTDIR", "./eig_resonator")
PALACE_BIN = os.environ["PALACE_BIN"]
NUM_CPUS   = int(os.environ.get("NSLOTS", "32"))
os.makedirs(OUTDIR, exist_ok=True)

TOTAL_LENGTH_MM  = float(os.environ.get("TOTAL_LENGTH_MM", "7.48345"))
# CPW cross-section. Default 10/6 um. Scaling both up (e.g. 20/12.25 um, keeps
# ~50 Ohm since w/(w+2g) barely moves) makes every feature ~2x larger and hence
# far easier for GMSH to resolve, at essentially unchanged eps_eff / frequency.
CPW_WIDTH_UM     = float(os.environ.get("CPW_WIDTH_UM", "10"))
CPW_GAP_UM       = float(os.environ.get("CPW_GAP_UM", "6"))
# Termination gap at the OPEN end. The qiskit-metal default is cpw_gap (6 um). A
# 6 um vacuum sliver in the ground plane is smaller than the background mesh
# target, so GMSH collapses it and Palace ties the resonator's "open" end to
# ground -> the line rings as a lambda/2 SHORT-SHORT mode (~8 GHz for 7.48 mm)
# instead of lambda/4 (~4 GHz). Widening the gap to ~20 um and fine-meshing the
# termination region (see below) keeps it a genuine open. 20 um still gives a
# sub-fF end capacitance, so the mode frequency is essentially unchanged.
OPEN_TERM_GAP_UM = float(os.environ.get("OPEN_TERM_GAP_UM", str(max(20.0, 3.0 * CPW_GAP_UM))))
# shift-and-invert target: must sit CLOSE to (ideally just below) the mode, or the
# shifted operator (K - sigma^2 M) is too indefinite and the multigrid-preconditioned
# FGMRES inside each eigen iteration diverges. lambda/4 (CPW on Si, eps_eff~6.35) is
# ~3.95 GHz for a 7.48 mm trace; open-end C and meander coupling pull it a bit lower.
STARTING_FREQ_HZ = float(os.environ.get("FSTART", "3.0e9"))
NUMBER_OF_FREQS  = int(os.environ.get("NFREQS", "3"))
SOLVER_ORDER     = int(os.environ.get("SOLVER_ORDER", "1"))   # order 2 -> ~5x the DoF (2M vs 320k) and less robust; bump later if Q needs it
RES_MIN_UM       = float(os.environ.get("RES_MIN_UM", "3.0"))
BG_MAX_UM        = 200.0

t0 = time.time()
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
import matplotlib
matplotlib.use("Agg", force=True)

design = designs.DesignPlanar({}, overwrite_enabled=True)
# chip sized tight around the bare resonator (no feedline -> no need for a wide chip;
# a smaller domain = fewer elements and pushes box/cavity modes above the resonator)
design.chips.main.size.size_x = '3.0mm'
design.chips.main.size.size_y = '2.4mm'
design.chips.main.size.size_z = '500um'
design.chips.main.size.center_x = '-0.2mm'
design.chips.main.size.center_y = '-0.65mm'
design.variables['cpw_width'] = f'{CPW_WIDTH_UM} um'
design.variables['cpw_gap'] = f'{CPW_GAP_UM} um'

OpenToGround(design, 'otg1', options=dict(chip='main', pos_x='-0.2mm', pos_y='-40um', orientation=180,
                                          width=f'{CPW_WIDTH_UM}um', gap=f'{CPW_GAP_UM}um',
                                          termination_gap=f'{OPEN_TERM_GAP_UM}um'))
ShortToGround(design, 'stg2', options=dict(chip='main', pos_x='0mm', pos_y='-1.25mm', orientation=0,
                                           width=f'{CPW_WIDTH_UM}um', gap=f'{CPW_GAP_UM}um'))
RouteMeander(design, 'resonator1', Dict(
    trace_width=f'{CPW_WIDTH_UM}um', trace_gap=f'{CPW_GAP_UM}um', total_length=f'{TOTAL_LENGTH_MM}mm',
    hfss_wire_bonds=False, fillet='99.9 um', lead=dict(start_straight='10um'),
    meander=Dict(spacing="200um", assymetry="1mm"),
    pin_inputs=Dict(start_pin=Dict(component='otg1', pin='open'),
                    end_pin=Dict(component='stg2', pin='short'))))
design.rebuild()
print("components:", list(design.components.keys()), f"  ({time.time()-t0:.1f}s)", flush=True)

# analytic guess for a quarter-wave CPW resonator on silicon (eps_eff ~ 6.45)
v = 3e8 / np.sqrt(6.45)
print(f"lambda/4 guess: f ~ {v / (4 * TOTAL_LENGTH_MM * 1e-3) / 1e9:.3f} GHz", flush=True)

from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.Utilities.Materials import MaterialInterface

eig = PALACE_Eigenmode_Simulation(
    name="eig_resonator", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "mesh_refinement": 0,
        "dielectric_material": "silicon",
        "starting_freq": STARTING_FREQ_HZ,
        "number_of_freqs": NUMBER_OF_FREQS,
        "solns_to_save": NUMBER_OF_FREQS,   # save every mode's field -> one ErealMag PNG per mode
        "solver_order": SOLVER_ORDER,
        "solver_tol": 1e-8,
        "solver_maxits": 500,
        "fillet_resolution": 12,
        "palace_dir": PALACE_BIN,
        "num_cpus": NUM_CPUS,
    })
eig.add_metallic(1)
eig.add_ground_plane()
# no ports: an isolated eigenmode of the bare resonator
# fine-mesh the trace AND both termination components: the open/short gaps live
# in the ground plane, not on `resonator1`, so without otg1/stg2 here the
# terminations get the 200 um background mesh and the open collapses to a short.
# taper distances are FIXED (not scaled with RES_MIN_UM) so a mesh-convergence
# sweep varies only the element size, not how fast the refinement falls off.
eig.fine_mesh_components(["resonator1", "otg1", "stg2"], min_size=RES_MIN_UM * 1e-6, max_size=BG_MAX_UM * 1e-6,
                         taper_dist_min=float(os.environ.get("TAPER_MIN_UM", "30")) * 1e-6,
                         taper_dist_max=float(os.environ.get("TAPER_MAX_UM", "120")) * 1e-6,
                         metals_only=False)
# belt-and-braces: a tight mesh box straddling the open termination gap
_box_min = max(1.5e-6, 0.25 * CPW_GAP_UM * 1e-6)
eig.fine_mesh_in_rectangle(-0.30e-3, -0.090e-3, -0.13e-3, -0.005e-3,
                           min_size=_box_min, max_size=BG_MAX_UM * 1e-6,
                           taper_dist_min=20e-6, taper_dist_max=80e-6)
eig.setup_EPR_interfaces(
    metal_air=MaterialInterface("Aluminium-Vacuum"),
    substrate_air=MaterialInterface("Silicon-Vacuum"),
    substrate_metal=MaterialInterface("Silicon-Aluminium"))

tp = time.time()
eig.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

ts = time.time()
eig.run()
print(f"SOLVE: {time.time()-ts:.1f} s", flush=True)

df = pd.read_csv(eig._output_data_dir + "/eig.csv")
df.columns = df.columns.str.strip()
fcol = [c for c in df.columns if c.startswith("Re{f}")][0]
qcol = [c for c in df.columns if c.startswith("Q")][0]
print("\n=== modes ===", flush=True)
for i, r in df.iterrows():
    print(f"  mode {i}: f = {r[fcol]:.5f} GHz,  Q = {r[qcol]:.3e}", flush=True)
df.to_csv(f"{OUTDIR}/modes.csv", index=False)
np.savez(f"{OUTDIR}/eig_resonator.npz", f_ghz=df[fcol].to_numpy(), Q=df[qcol].to_numpy())
print(f"\nTOTAL wall: {time.time()-t0:.1f} s", flush=True)
