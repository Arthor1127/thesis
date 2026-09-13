"""Driven S21 validation of what create_port_2_conds's impedance_C actually
realizes -- see shuntc_validation_design.py's docstring for the full
motivation and the exact analytic S21(f) formula for an ideal shunt
capacitor. No resonator, no mode-finding: a plain broadband driven sweep,
comparing simulated S21 against the closed-form prediction at every
frequency, and inverting to get the effective realized C(f).

Env: CC_FF (target capacitance in femtofarad), PALACE_BIN, [NSLOTS=16],
     [OUTDIR=./out_shuntc], [TAG], [FSTART_GHZ=1], [FEND_GHZ=15], [NPTS=141]
Writes: {OUTDIR}/results/{TAG}.json, {OUTDIR}/results/{TAG}_s21.npz
Cluster job only.
"""
import os
import sys
import json
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

CC_FF        = float(os.environ["CC_FF"])
OUTDIR       = os.environ.get("OUTDIR", "./out_shuntc")
PALACE_BIN   = os.environ["PALACE_BIN"]
NCPU         = int(os.environ.get("NSLOTS", "16"))
TAG          = os.environ.get("TAG", f"shuntc_{CC_FF:g}fF")
FSTART_GHZ   = float(os.environ.get("FSTART_GHZ", "1.0"))
FEND_GHZ     = float(os.environ.get("FEND_GHZ", "15.0"))
N_POINTS     = int(os.environ.get("NPTS", "141"))
os.makedirs(f"{OUTDIR}/results", exist_ok=True)

Cc_F = CC_FF * 1e-15
t0 = time.time()
from qiskit_metal import designs  # noqa: F401
import matplotlib
matplotlib.use("Agg", force=True)

import shuntc_validation_design as scd

print(f"[{TAG}] Cc={CC_FF}fF ncpu={NCPU} sweep {FSTART_GHZ}-{FEND_GHZ}GHz x{N_POINTS} "
     f"(no resonance -- plain broadband sweep, no PROM needed)", flush=True)

design = scd.build_design(Cc_F)
print(f"geometry built ({time.time()-t0:.1f}s)  chip={design.chips.main.size.size_x} "
      f"x {design.chips.main.size.size_y}", flush=True)

from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation

driven = PALACE_Driven_Simulation(
    name=f"shuntc_{TAG}", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    user_options={
        "dielectric_material": "silicon", "solns_to_save": 0,
        "solver_order": 1, "solver_tol": 1e-6, "solver_maxits": 500,
        "solver_initial_guess": True, "fillet_resolution": 12,
        "palace_dir": PALACE_BIN, "num_cpus": NCPU,
    })
driven.add_metallic(1)
driven.add_ground_plane()
driven.create_port_CPW_on_Route("TL1", pin_name="start", len_launch=30e-6, impedance_R=50)
driven.create_port_CPW_on_Route("TL2", pin_name="end", len_launch=30e-6, impedance_R=50)
# PORT_R_OHM: diagnostic knob added after the first pass showed the
# "realized C" was ~600fF regardless of the requested Cc (2.88-100fF) and
# NOT frequency-independent (35% variation across the sweep) -- i.e. NOT
# behaving like a lumped capacitor at all, even though config.json faithfully
# recorded the correct C value. Hypothesis: R=0 doesn't mean "pure reactance,
# zero series resistance" to Palace's Robin lumped-port solver -- it may hit
# a degenerate/different code path (e.g. treated as a hard short) instead.
# Testing a tiny-but-nonzero R instead of exactly 0.
PORT_R_OHM = float(os.environ.get("PORT_R_OHM", "0"))
driven.create_port_2_conds("JT", "second_end", "GND1", "short",
                           rect_width=20e-6, impedance_R=PORT_R_OHM, impedance_C=Cc_F)
driven.set_port_excitation(1)

lateral_min_um = 0.5 * min(scd.CPW_WIDTH_UM, scd.CPW_GAP_UM, scd.COUPLING_GAP_UM)
driven.fine_mesh_components(["TL1", "TL2", "JT", "OTG1", "OTG2", "GND1"],
                            min_size=lateral_min_um * 1e-6, max_size=150e-6,
                            taper_dist_min=16e-6, taper_dist_max=40e-6, metals_only=False)
print(f"mesh: lateral_min={lateral_min_um:.3f}um", flush=True)

fs, fe = FSTART_GHZ * 1e9, FEND_GHZ * 1e9
driven.set_freq_values(freq_start=fs, freq_end=fe, freq_step=(fe - fs) / max(N_POINTS - 1, 1))

tp = time.time()
driven.prepare_simulation()
print(f"MESH+PREP: {time.time()-tp:.1f} s", flush=True)

# No resonance to resolve -- a plain, coarse-relative-to-run_qc_sim adaptive
# sweep (or even a discrete one) is enough; PROM defaults are fine here.
cfg_path = driven._sim_config
with open(cfg_path) as f:
    cfg = json.load(f)
cfg["Solver"]["Driven"]["AdaptiveTol"] = 1e-4
cfg["Solver"]["Driven"]["AdaptiveMaxSamples"] = 30
with open(cfg_path, "w") as f:
    json.dump(cfg, f, indent=2)

ts = time.time()
driven.run()
solve_s = time.time() - ts
print(f"SOLVE: {solve_s:.1f} s", flush=True)

s = driven.retrieve_data()
if s is None:
    raise SystemExit("retrieve_data() returned None -- the Palace solve failed.")

import numpy as np
freqs = np.asarray(s["freqs"]); S21 = np.asarray(s["S21"])
np.savez(f"{OUTDIR}/results/{TAG}_s21.npz", freq_hz=freqs, S21=S21, Cc_F=Cc_F)

S21_pred = scd.S21_analytic(freqs, Cc_F)
C_realized = scd.C_from_S21(freqs, S21)

# Use the mid-band point (least sensitive to any residual mesh/discretization
# noise near the sweep edges) as the headline "realized C" estimate, but save
# the full frequency-resolved array too -- if realized C is frequency-
# INDEPENDENT, that confirms it behaves as a genuine lumped capacitor (just
# possibly scaled); if it drifts with frequency, something else is going on
# (e.g. a real parasitic inductance/resistance contaminating the port).
mid = len(freqs) // 2
C_realized_mid_fF = float(np.real(C_realized[mid])) * 1e15
C_ratio_mid = float(np.real(C_realized[mid])) / Cc_F

row = {
    "tag": TAG, "Cc_target_fF": CC_FF, "Z0_ohm": scd.Z0_OHM, "port_R_ohm": PORT_R_OHM,
    "C_realized_mid_fF": C_realized_mid_fF, "C_ratio_mid": C_ratio_mid,
    "freq_mid_ghz": float(freqs[mid] / 1e9),
    "S21_mid_sim": float(np.real(S21[mid])), "S21_mid_sim_imag": float(np.imag(S21[mid])),
    "S21_mid_pred": float(np.real(S21_pred[mid])), "S21_mid_pred_imag": float(np.imag(S21_pred[mid])),
    "C_realized_std_over_mean_pct": float(100 * np.std(np.real(C_realized)) / np.mean(np.real(C_realized))),
    "fstart_ghz": FSTART_GHZ, "fend_ghz": FEND_GHZ, "n_points": N_POINTS,
    "solve_s": solve_s, "wall_s": time.time() - t0,
}
with open(f"{OUTDIR}/results/{TAG}.json", "w") as f:
    json.dump(row, f, indent=2)

print(json.dumps(row, indent=2), flush=True)
print(f"\nC_realized(f) [fF]: {(C_realized*1e15).tolist()}", flush=True)
print(f"TOTAL wall: {time.time()-t0:.1f} s", flush=True)
