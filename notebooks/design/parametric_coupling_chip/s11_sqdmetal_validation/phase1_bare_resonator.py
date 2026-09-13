"""Phase 1 of the PalaceForCQED S11 cross-check (see HANDOFF.md): feedline +
ONE bare hanger resonator (no qubit yet), dimensioned to match resonator
"element 0" of PalaceForCQED's own reference geometry
(palace_cqed_reference/2DQv8_eb4_edit.geo), reverse-engineered directly from
that .geo file's raw point coordinates (see scratchpad parse/plot scripts --
not re-derivable from this script alone, so the numbers below are asserted,
not re-computed):

  - CPW cross-section: trace_width=20um, trace_gap=10um -- same on the
    feedline and on all 4 hanger resonators in the reference (self-consistent
    single design rule, confirmed by x-crossing profiling of both the
    feedline (elements 4/5) and every resonator (elements 0-3)).
  - Resonator total length ~4120um (shoelace area / (2*trace_gap) of element
    0's polygon) -- a quarter-wave (one end open/coupled, one end shorted) at
    this length on Si (eps_r=11.49) lands at ~7.2GHz, matching the observed
    7.1-7.7GHz dip cluster in driven_new.svg (the two lower-frequency driven
    windows in run.jl, 3.83-3.93/4.35-4.45GHz, are the qubit modes, not more
    resonators -- confirmed since all 4 resonators came back the same length).
  - Coupling gap to feedline ~115um at closest approach (element 0's nearest
    edge to the feedline) -- a real, fairly weak/loose coupling, consistent
    with driven_new.svg's dips being shallow (-6 to -9dB) rather than deep.
  - The reference's real fold pattern (legs 20um apart within a U-turn, 40um
    between different U-turns) is NOT reproduced exactly -- this uses
    qiskit_metal's own snake/meander routing (same pattern as
    readout_resonators.py) to hit the right TOTAL LENGTH and coupling gap,
    which set the resonant frequency and Qc; the exact fold geometry is a
    second-order effect on self-coupling, not attempted here.
  - Quarter-wave assumption (shorted far end): not directly confirmed from
    the .geo file's end connectivity (that would need tracing which polygon,
    if any, the meander's far tip merges into) -- inferred from the length-vs-
    observed-frequency match. Flagged, not verified at the byte level.

Single EXCITATION port only (matching PalaceForCQED's own config: ONE lumped
port, R=50 Ohm, at one end of the feedline; the other end is a true open
circuit) -- this is a 1-port |S11| device, not a 2-port S21 line.

This script only builds the geometry and (if MESH=1) runs prepare_simulation()
(GMSH meshing -- NOT a PALACE solve, safe to run locally per this project's
"local-machine-heavy-fem-jobs-unsafe" rule, which explicitly allows local
mesh-size inspection). It does NOT call driven.run() -- no cluster job is
submitted by this script. See s11_sqdmetal_validation/HANDOFF.md.

Env: [MESH=0], [OUTDIR=./out_phase1]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # parametric_coupling_chip/, for readout_resonators

import matplotlib
matplotlib.use("Agg", force=True)

from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
from qiskit_metal.qlibrary.couplers.coupled_line_tee import CoupledLineTee

from readout_resonators import build_snake_anchors, resonator_overlap_um2

# ---- dimensions measured directly from palace_cqed_reference/2DQv8_eb4_edit.geo ----
CPW_WIDTH_UM = 20.0
CPW_GAP_UM = 10.0
TARGET_RES_LENGTH_MM = 4.12          # element 0's shoelace-derived length
COUPLING_SPACE_UM = 115.0            # element 0's closest approach to the feedline
COUPLING_LENGTH_UM = 400.0           # element 0's near leg length (parallel run)
COUPLING_BUDGET_UM = 700.0           # coupling_length_um + down_length_um (readout_resonators.py convention)
FILLET_RADIUS_UM = 50.0
SNAKE_FILLET_UM = 100.0

FEEDLINE_HALF_LEN_UM = 2000.0        # ~matches the reference's ~3.9mm main feedline run
# CHIP_SIZE_X_MM: trimmed margin (900->300um) plus a smaller x_prop air-pad
# below -- the feedline direction had far more boundary clearance than it
# needed. CHIP_SIZE_Y_MM: raised 3.0->4.0mm -- the resonator's snake reaches
# to about y=-1.3mm and needs more room to the y boundary than it had.
CHIP_SIZE_X_MM = 2.5 * (FEEDLINE_HALF_LEN_UM + COUPLING_SPACE_UM + 300.0) / 1000.0
CHIP_SIZE_Y_MM = 3.0

down_length_um = COUPLING_BUDGET_UM - COUPLING_LENGTH_UM
assert down_length_um > 0

design = designs.DesignPlanar({}, overwrite_enabled=True)
design.chips.main.size.size_x = f'{CHIP_SIZE_X_MM}mm'
design.chips.main.size.size_y = f'{CHIP_SIZE_Y_MM}mm'
design.chips.main.size.size_z = '500um'
design.chips.main.size.center_x = '0mm'
design.chips.main.size.center_y = '-0.5mm'
design.variables['cpw_width'] = f'{CPW_WIDTH_UM}um'
design.variables['cpw_gap'] = f'{CPW_GAP_UM}um'

# Feedline runs along x (matching readout_resonators.py's proven LP1-TL1-CLT-
# TL2-LP2 layout exactly -- a y-axis feedline hit qiskit_metal's known
# RouteAnchors/RoutePathfinder obstruction-check fragility (see
# readout_resonators.build_design's docstring on the same 'NoneType' object
# has no attribute 'coords' error); absolute orientation has no physical
# meaning, so just reuse the axis this codebase's tooling is proven on).
# One true open end (OTG_LEFT), one end that becomes a lumped excitation port
# (OTG_RIGHT is just a pin anchor for RoutePathfinder -- create_port_CPW_on_Route
# overrides it with the real port, same convention as check_line_impedance.py's
# line_only mode).
OpenToGround(design, 'OTG_LEFT', options=dict(
    chip='main', pos_x=f'-{FEEDLINE_HALF_LEN_UM}um', pos_y='0um', orientation='180',
    width='cpw_width', gap='cpw_gap'))
OpenToGround(design, 'OTG_RIGHT', options=dict(
    chip='main', pos_x=f'{FEEDLINE_HALF_LEN_UM}um', pos_y='0um', orientation='0',
    width='cpw_width', gap='cpw_gap'))

CoupledLineTee(design, 'CLT', options=dict(
    pos_x='0um', pos_y='0um', orientation='0', chip='main',
    prime_width='cpw_width', prime_gap='cpw_gap',
    second_width='cpw_width', second_gap='cpw_gap',
    coupling_space=f'{COUPLING_SPACE_UM}um', coupling_length=f'{COUPLING_LENGTH_UM}um',
    down_length=f'{down_length_um}um', fillet=f'{FILLET_RADIUS_UM}um', mirror=False,
    open_termination=True))

RoutePathfinder(design, 'TL1', options=dict(
    chip='main', trace_width='cpw_width', trace_gap='cpw_gap', fillet=f'{FILLET_RADIUS_UM}um',
    lead=dict(end_straight='0.1mm'),
    pin_inputs=Dict(start_pin=Dict(component='OTG_LEFT', pin='open'),
                    end_pin=Dict(component='CLT', pin='prime_start'))))
RoutePathfinder(design, 'TL2', options=dict(
    chip='main', trace_width='cpw_width', trace_gap='cpw_gap', fillet=f'{FILLET_RADIUS_UM}um',
    lead=dict(start_straight='0.1mm'),
    pin_inputs=Dict(start_pin=Dict(component='CLT', pin='prime_end'),
                    end_pin=Dict(component='OTG_RIGHT', pin='open'))))

design.rebuild()

# Snake the resonator off CLT's second_end to hit TARGET_RES_LENGTH_MM total,
# same boustrophedon helper readout_resonators.py uses for its own sweep.
consumed_mm = (COUPLING_LENGTH_UM + down_length_um) * 1e-3
remaining_mm = TARGET_RES_LENGTH_MM - consumed_mm
assert remaining_mm > 0.2, f"remaining_mm={remaining_mm:.3f} too short"

se = design.components['CLT'].pins['second_end']['middle']
pts_um, heading = build_snake_anchors(se, remaining_mm, x_max_um=600.0)
orientation = {'EAST': '0', 'WEST': '180', 'DOWN': '270'}[heading]
end_x_um, end_y_um = pts_um[-1]

ShortToGround(design, 'res_short', options=dict(
    pos_x=f'{end_x_um}um', pos_y=f'{end_y_um}um', orientation=orientation,
    chip='main', width='cpw_width'))

from collections import OrderedDict
import numpy as np
anchors = OrderedDict()
for i, (x_um, y_um) in enumerate(pts_um[:-1]):
    anchors[i] = np.array([x_um, y_um]) * 1e-3

RouteAnchors(design, 'readout_1', options=dict(
    pin_inputs=dict(
        start_pin=dict(component='CLT', pin='second_end'),
        end_pin=dict(component='res_short', pin='short'),
    ),
    fillet=f'{SNAKE_FILLET_UM}um',
    anchors=anchors,
    lead=dict(start_straight='2um', end_straight='2um'),
    trace_width='cpw_width', trace_gap='cpw_gap',
))

design.rebuild()
overlap = resonator_overlap_um2(design, other_names=("TL1", "TL2", "CLT"))
print(f"geometry built. resonator/feedline overlap = {overlap:.4f} um^2 (should be 0)")
print(f"chip size: {design.chips.main.size.size_x} x {design.chips.main.size.size_y}")
print(f"resonator snake heading={heading}, {len(anchors)} anchors, target_length={TARGET_RES_LENGTH_MM}mm "
      f"(coupling {COUPLING_LENGTH_UM}um + down {down_length_um}um + snake {remaining_mm*1e3:.0f}um)")

OUTDIR = os.environ.get("OUTDIR", "./out_phase1")
os.makedirs(OUTDIR, exist_ok=True)

# ---- 2D layout preview (qgeometry polygons, colored by component) ----
import shapely
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(10, 14))
pt = design.qgeometry.tables["path"]
pp = design.qgeometry.tables["poly"]
comp_names = {v: k for k, v in {n: c.id for n, c in design.components.items()}.items()}
cmap = plt.cm.tab10
color_for = {}
def get_color(cid):
    if cid not in color_for:
        color_for[cid] = cmap(len(color_for) % 10)
    return color_for[cid]

for _, r in pt.iterrows():
    w = float(r.get("width") or 0.0) or 0.01
    poly = r["geometry"].buffer(w / 2.0, cap_style=2, join_style=2)
    xs, ys = poly.exterior.xy
    ax.fill(xs, ys, color=get_color(r["component"]), alpha=0.8,
            label=comp_names.get(r["component"], str(r["component"])))
for _, r in pp.iterrows():
    geom = r["geometry"]
    if geom.geom_type == "Polygon":
        xs, ys = geom.exterior.xy
        ax.fill(xs, ys, color=get_color(r["component"]), alpha=0.5)

handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), loc='upper right', fontsize=8)
ax.set_aspect('equal')
ax.set_title(f"Phase 1: feedline + bare resonator (matches element 0 dims)\n"
             f"CPW {CPW_WIDTH_UM}/{CPW_GAP_UM}um, L_res={TARGET_RES_LENGTH_MM}mm, "
             f"coupling_space={COUPLING_SPACE_UM}um, coupling_length={COUPLING_LENGTH_UM}um")
ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
fig.savefig(f"{OUTDIR}/layout_preview.png", dpi=150, bbox_inches='tight')
print(f"saved {OUTDIR}/layout_preview.png")

fig2, ax2 = plt.subplots(figsize=(8, 6))
for _, r in pt.iterrows():
    w = float(r.get("width") or 0.0) or 0.01
    poly = r["geometry"].buffer(w / 2.0, cap_style=2, join_style=2)
    xs, ys = poly.exterior.xy
    ax2.fill(xs, ys, color=get_color(r["component"]), alpha=0.8,
             label=comp_names.get(r["component"], str(r["component"])))
handles2, labels2 = ax2.get_legend_handles_labels()
by_label2 = dict(zip(labels2, handles2))
ax2.legend(by_label2.values(), by_label2.keys(), loc='upper right', fontsize=8)
ax2.set_xlim(-0.3, 0.5)
ax2.set_ylim(-0.6, 0.15)
ax2.set_aspect('equal')
ax2.set_title(f"Zoom: CoupledLineTee coupling region (coupling_space={COUPLING_SPACE_UM}um)")
ax2.set_xlabel("x (mm)"); ax2.set_ylabel("y (mm)")
fig2.savefig(f"{OUTDIR}/coupling_zoom.png", dpi=200, bbox_inches='tight')
print(f"saved {OUTDIR}/coupling_zoom.png")

if os.environ.get("MESH", "0") == "1":
    from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation
    driven = PALACE_Driven_Simulation(
        name="phase1_bare_res", metal_design=design, sim_parent_directory=OUTDIR,
        mode="simPC", meshing="GMSH", create_files=True,
        user_options={
            "dielectric_material": "silicon", "solns_to_save": 0,
            "solver_order": 4, "solver_tol": 1e-8, "solver_maxits": 300,
            "solver_initial_guess": True, "fillet_resolution": 12,
            "palace_dir": "/dev/null", "num_cpus": 1,
        })
    driven.add_metallic(1)
    driven.add_ground_plane()
    driven.set_farfield("absorbing")
    # Simple in-line port (no LaunchpadWirebond taper, no interdigitated
    # coupling capacitor) -- the reference's real excitation port has both
    # (see palace_cqed_reference elements 6/10/12/14/16/17/18/19/20/21). Kept
    # simple on the basis of this repo's own earlier check_line_impedance.py
    # finding that a launchpad taper's contribution to a bare line's behavior
    # was small -- NOT independently re-verified on this exact geometry, so
    # treat as a carried-over assumption, not a confirmed result here.
    driven.create_port_CPW_on_Route("TL2", pin_name="end", len_launch=CPW_GAP_UM * 1e-6, impedance_R=50.0)
    driven.set_port_excitation(1)

    # Coarser base mesh on TL/CLT/readout_1 (min_size factor 0.5->0.75, i.e.
    # 5um->20um) with a wider, gentler taper out to the background max_size
    # (also raised 100->200um) -- relying on AMR (below) to refine further
    # only where the solution actually needs it, instead of a uniformly fine
    # hand-set base mesh everywhere along the trace.
    MIN_SIZE_FACTOR = float(os.environ.get("MIN_SIZE_FACTOR", "1.2"))
    MAX_SIZE_UM = float(os.environ.get("MAX_SIZE_UM", "190.0"))
    TAPER_DIST_MIN_UM = float(os.environ.get("TAPER_DIST_MIN_UM", "20.0"))
    TAPER_DIST_MAX_UM = float(os.environ.get("TAPER_DIST_MAX_UM", "400.0"))
    lateral_min_um = MIN_SIZE_FACTOR * min(CPW_WIDTH_UM, CPW_GAP_UM)
    driven.fine_mesh_components(["TL1", "TL2", "CLT", "readout_1"],
                                min_size=lateral_min_um * 1e-6, max_size=MAX_SIZE_UM * 1e-6,
                                taper_dist_min=TAPER_DIST_MIN_UM * 1e-6, taper_dist_max=TAPER_DIST_MAX_UM * 1e-6,
                                metals_only=False)

    AMR_ITERS = int(os.environ.get("AMR_ITERS", "4"))
    if AMR_ITERS > 0:
        AMR_TOL = float(os.environ.get("AMR_TOL", "1e-2"))
        AMR_MAX_DOFS = int(os.environ.get("AMR_MAX_DOFS", "3000000"))
        driven.enable_mesh_refinement(num_iterations=AMR_ITERS, max_DoFs=AMR_MAX_DOFS,
                                      tolerance=AMR_TOL, Dorfler_marking_fraction=0.7,
                                      save_iterations_data=True, save_iterations_mesh=False)

    driven.set_zBoundary_as_proportion(1.0, 0.3)
    # x_prop trimmed (0.1->0.05): the feedline had more boundary clearance
    # than needed. y_prop raised (0.1->0.15) on top of CHIP_SIZE_Y_MM's own
    # increase above: more room between the resonator's snake and the
    # boundary. x and y are independent knobs here since the feedline
    # extends along x and the resonator's snake extends along y.
    driven.set_xBoundary_as_proportion(0.05)
    driven.set_yBoundary_as_proportion(0.15)
    driven.set_freq_values(freq_start=7.0e9, freq_end=7.8e9, freq_step=0.1e9)
    driven.prepare_simulation()
    print(f"mesh prepared in {OUTDIR} (GMSH only, no PALACE solve run); "
          f"base lateral_min={lateral_min_um}um, AMR_ITERS={AMR_ITERS}")
