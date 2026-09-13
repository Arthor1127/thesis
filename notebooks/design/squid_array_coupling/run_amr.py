"""
AMR run for the SQUID-array + readout-resonator design, then a per-junction
EPR extraction for the readout-resonator mode.

The SQUID array is a ~1.1mm slow-wave line whose standing-wave "phonon" modes
follow the acoustic dispersion (C = island-to-ground capacitance, which is
large here - the loop metal sits a few um from the ground plane):

    omega_j = 2 sin(k_j/2) / sqrt( L_J (4 C_J sin^2(k_j/2) + C_g) ),  k_j = j*pi/(N+1)

Band top ~ 1/sqrt(L_J), so a large L_J (>=400 nH) drops the whole band below
the ~8 GHz readout resonator. The one mode whose field lives on the CPW meander
- the resonator - is picked by identify_resonator_mode() (sum_j |p_j| << 1),
NOT by frequency. Its full per-junction participation vector is written to
epr_result.json for extract_c_array.py.

Env overrides (all optional):
    PALACE_BIN        (required)  path to the palace binary
    NSLOTS            (16)        MPI ranks
    N_SQUIDS_OVR      (build_design.N_SQUIDS)  loops in the array
    L_J_NH           (build_design.L_J_NH)    inductance of every junction, nH
    AMR_ITERS         (5)         AMR iterations
    AMR_MAX_DOFS      (15000000)  stop refining past this many DoF (RAM guard)
    AMR_DORFLER       (0.6)       Dorfler marking fraction
    STARTING_FREQ_HZ  (6.5e9)     eigen-search start (above the phonon band)
    NUMBER_OF_FREQS   (3)         modes to solve
    RESONATOR_MIN_UM  (4.0)       readout_1 base fine-mesh min size
    ARRAY_MIN_UM      (3.5)       array_1  base fine-mesh min size
    RESONATOR_TAPER_UM(40.0)      readout_1 fine-zone taper distance
    ARRAY_TAPER_UM    (15.0)      array_1  fine-zone taper distance
Args:
    sys.argv[1]  (default runs/amr_run1)  output directory
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_design as bd
from build_design import (build_design, identify_resonator_mode,
                          CJ_J_FF, RES_LENGTH_MM, JJ_PER_LOOP, ARRAY_COMPONENT_NAME)


def _env(name, default, cast=float):
    v = os.environ.get(name)
    return default if v is None or v == "" else cast(v)


OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "runs/amr_run1"
PALACE_BIN = os.environ["PALACE_BIN"]
NUM_CPUS = _env("NSLOTS", 16, int)
N_SQUIDS = _env("N_SQUIDS_OVR", bd.N_SQUIDS, int)
L_J_nH = _env("L_J_NH", bd.L_J_NH)
AMR_ITERS = _env("AMR_ITERS", 5, int)
AMR_MAX_DOFS = _env("AMR_MAX_DOFS", 15_000_000, int)
AMR_DORFLER = _env("AMR_DORFLER", 0.6)
STARTING_FREQ = _env("STARTING_FREQ_HZ", 6.5e9)
NUMBER_OF_FREQS = _env("NUMBER_OF_FREQS", 3, int)
RES_MIN_UM = _env("RESONATOR_MIN_UM", 4.0)
ARR_MIN_UM = _env("ARRAY_MIN_UM", 3.5)
RES_TAPER_UM = _env("RESONATOR_TAPER_UM", 40.0)
ARR_TAPER_UM = _env("ARRAY_TAPER_UM", 15.0)
os.makedirs(OUTDIR, exist_ok=True)

N_PORTS = JJ_PER_LOOP * N_SQUIDS

import time
t0 = time.time()
print(f"[1/6] Building design (n_squids={N_SQUIDS}, L_j={L_J_nH}nH, C_j={CJ_J_FF}fF)...", flush=True)
design = build_design(n_squids=N_SQUIDS, L_j_nH=L_J_nH)
design.rebuild()
n_jj = len(design.qgeometry.tables['junction'].loc[
    design.qgeometry.tables['junction']['component'] == design.components[ARRAY_COMPONENT_NAME].id])
assert n_jj == N_PORTS, f"expected {N_PORTS} junction rows, got {n_jj}"
print(f"      OK ({time.time()-t0:.1f}s). Components: {list(design.components.keys())}, "
      f"{n_jj} junctions", flush=True)

print("[2/6] Importing SQDMetal PALACE eigenmode module...", flush=True)
from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.Utilities.Materials import MaterialInterface
import matplotlib
matplotlib.use("Agg", force=True)
print("      OK", flush=True)

print(f"[3/6] Eigenmode sim + {N_PORTS} junction ports "
      f"(solve {NUMBER_OF_FREQS} modes from {STARTING_FREQ/1e9:.2f} GHz)...", flush=True)
t1 = time.time()
user_defined_options = {
    "mesh_refinement": 0,
    "dielectric_material": "silicon",
    "starting_freq": STARTING_FREQ,
    "number_of_freqs": NUMBER_OF_FREQS,
    "solns_to_save": NUMBER_OF_FREQS,
    "solver_order": 1,
    "solver_tol": 1.0e-5,
    "solver_maxits": 300,
    "fillet_resolution": 10,
    "palace_dir": PALACE_BIN,
    "num_cpus": NUM_CPUS,
}
eigen_sim = PALACE_Eigenmode_Simulation(
    name="amr_array", metal_design=design, sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", user_options=user_defined_options, create_files=True)
eigen_sim.add_metallic(1)
eigen_sim.add_ground_plane()
for k in range(N_PORTS):
    eigen_sim.create_port_JosephsonJunction(
        ARRAY_COMPONENT_NAME, junction_index=k,
        L_J=L_J_nH * 1e-9, C_J=CJ_J_FF * 1e-15)
print(f"      OK ({time.time()-t1:.1f}s). Ports: {len(eigen_sim._ports)}", flush=True)

print(f"[4/6] Base fine mesh: readout_1 min={RES_MIN_UM}um taper={RES_TAPER_UM}um, "
      f"{ARRAY_COMPONENT_NAME} min={ARR_MIN_UM}um taper={ARR_TAPER_UM}um "
      "(user 2026-08-28: AMR was bisecting resonator + array repeatedly - "
      "seed them finer)...", flush=True)
t2 = time.time()
eigen_sim.fine_mesh_components(
    ["readout_1"], min_size=RES_MIN_UM * 1e-6, max_size=150e-6,
    taper_dist_min=RES_TAPER_UM * 1e-6, metals_only=False)
eigen_sim.fine_mesh_components(
    [ARRAY_COMPONENT_NAME], min_size=ARR_MIN_UM * 1e-6, max_size=150e-6,
    taper_dist_min=ARR_TAPER_UM * 1e-6, metals_only=False)
eigen_sim.setup_EPR_interfaces(
    metal_air=MaterialInterface("Aluminium-Vacuum"),
    substrate_air=MaterialInterface("Silicon-Vacuum"),
    substrate_metal=MaterialInterface("Silicon-Aluminium"))
print(f"      OK ({time.time()-t2:.1f}s)", flush=True)

print(f"[5/6] AMR: {AMR_ITERS} iters, tol=5e-3, Dorfler={AMR_DORFLER}, "
      f"max_DoFs={AMR_MAX_DOFS:,} (RAM guard), saving per-iteration data...", flush=True)
eigen_sim.enable_mesh_refinement(
    num_iterations=AMR_ITERS, max_DoFs=AMR_MAX_DOFS, tolerance=5e-3,
    Dorfler_marking_fraction=AMR_DORFLER,
    save_iterations_data=True, save_iterations_mesh=True, nonconformal=False)

print("[6/6] prepare_simulation() + run()...", flush=True)
t3 = time.time()
eigen_sim.prepare_simulation()
eigen_sim.run()
print(f"      OK ({time.time()-t3:.1f}s)", flush=True)

# ---------- per-iteration convergence of the resonator frequency ----------
import glob
import pandas as pd

OF = eigen_sim._output_data_dir
iter_dirs = sorted(glob.glob(os.path.join(OF, "iteration*")))
conv = []
for d in iter_dirs + [OF]:   # each AMR iteration, then the final postpro dir
    ep = os.path.join(d, "port-EPR.csv")
    ec = os.path.join(d, "eig.csv")
    if not (os.path.exists(ep) and os.path.exists(ec)):
        continue
    P = pd.read_csv(ep).to_numpy()[:, 1:]
    E = pd.read_csv(ec).to_numpy()
    fz = E[:, 1] * 1e9
    try:
        ri, _ = identify_resonator_mode(P, fz)
        conv.append({"where": os.path.basename(d) or "final",
                     "f_ghz": float(fz[ri] / 1e9),
                     "sum_abs_p": float(np.abs(P[ri]).sum()),
                     "n_dof_note": None})
    except RuntimeError as e:
        conv.append({"where": os.path.basename(d) or "final", "error": str(e)})

# dedupe (final dir often == last iteration dir)
seen_f = set()
conv_u = []
for c in conv:
    key = round(c.get("f_ghz", -1), 6)
    if key in seen_f:
        continue
    seen_f.add(key)
    conv_u.append(c)
for i in range(1, len(conv_u)):
    if "f_ghz" in conv_u[i] and "f_ghz" in conv_u[i - 1]:
        conv_u[i]["rel_change_vs_prev"] = abs(
            conv_u[i]["f_ghz"] - conv_u[i - 1]["f_ghz"]) / conv_u[i]["f_ghz"]

converged = (len(conv_u) >= 2 and "rel_change_vs_prev" in conv_u[-1]
             and conv_u[-1]["rel_change_vs_prev"] < 5e-3)

# ---------- final resonator-mode EPR ----------
port_epr = eigen_sim.retrieve_mode_port_EPR()
mesh_sizes = eigen_sim.retrieve_simulation_sizes()
mat = port_epr["mat_mode_port"]
fz = port_epr["eigenfrequencies"]
Q = port_epr["loaded_Q"]
res_idx, res_info = identify_resonator_mode(mat, fz)

p_signed = mat.sum(axis=1)
p_abs = np.abs(mat).sum(axis=1)

omega0 = 2.0 * math.sqrt(1.0 / (L_J_nH * 1e-9 * CJ_J_FF * 1e-15))

resonator = {
    "m": int(res_idx + 1),
    "f_ghz": float(fz[res_idx] / 1e9),
    "Q": float(Q[res_idx]),
    "p_array_signed": float(p_signed[res_idx]),
    "p_array_abs_sum": float(p_abs[res_idx]),
    "p_per_junction": [float(x) for x in mat[res_idx]],
}
if mat.shape[1] == N_PORTS:
    resonator["p_per_loop"] = [
        float(mat[res_idx, 2 * i] + mat[res_idx, 2 * i + 1]) for i in range(N_SQUIDS)]

result = {
    "status": "success",
    "n_squids": N_SQUIDS,
    "l_j_nh": L_J_nH,
    "cj_j_ff": CJ_J_FF,
    "res_length_mm": RES_LENGTH_MM,
    "n_inductive_ports": int(mat.shape[1]),
    "amr": {
        "iterations_requested": AMR_ITERS,
        "max_dofs": AMR_MAX_DOFS,
        "dorfler": AMR_DORFLER,
        "base_mesh_um": {"resonator_min": RES_MIN_UM, "array_min": ARR_MIN_UM},
        "convergence": conv_u,
        "frequency_converged_lt_0.5pct": bool(converged),
    },
    "final_mesh": {"elements": mesh_sizes.get("MeshElements"), **mesh_sizes},
    "resonator_id": res_info,
    "all_modes": [{
        "m": i + 1, "f_ghz": float(fz[i] / 1e9), "Q": float(Q[i]),
        "p_array_signed": float(p_signed[i]), "p_array_abs_sum": float(p_abs[i]),
        "kind": "resonator" if i == res_idx else "phonon/box",
    } for i in range(len(fz))],
    "resonator_mode": resonator,
    "phonon_dispersion": {
        "formula": "omega_j = 2*sqrt(1/(L_J*C)) * sin(j*pi/(2*(N+1)))  [acoustic, C=C_g]",
        "omega0_over_2pi_ghz_using_CJ": omega0 / (2 * math.pi) / 1e9,
    },
    "wall_time_s": time.time() - t0,
}
out_path = os.path.join(OUTDIR, "epr_result.json")
json.dump(result, open(out_path, "w"), indent=2)

print(f"\n=== AMR RUN COMPLETE in {time.time()-t0:.1f}s ===", flush=True)
print(f"Final mesh: {mesh_sizes}", flush=True)
print("Resonator-frequency convergence:", flush=True)
for c in conv_u:
    print(f"  {c.get('where'):>12}  f={c.get('f_ghz','?'):.5f} GHz  "
          f"d={c.get('rel_change_vs_prev','-')}", flush=True)
print(f"CONVERGED (<0.5%): {converged}", flush=True)
print(f"RESONATOR: m={res_idx+1}, f={fz[res_idx]/1e9:.4f} GHz, "
      f"sum|p_j|={p_abs[res_idx]:.3e}, signed={p_signed[res_idx]:.3e}", flush=True)
print(f"Result -> {out_path}", flush=True)
