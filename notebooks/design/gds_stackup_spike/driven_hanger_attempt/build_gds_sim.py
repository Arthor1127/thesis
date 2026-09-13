"""Feasibility test: feed PalaceForCQED's real .geo elements (resonator +
feedline, converted 1:1 to GDS by geo_to_gds.py) into SQDMetal via
gds_design=, using PALACE_Driven_Simulation directly (gds-stackup-pipeline's
own build.py only wraps eigenmode so far -- extending its pattern here by
hand for driven). Geometry + mesh only, no run(), no cluster job."""
import os
import matplotlib
matplotlib.use("Agg", force=True)

GDS_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "phase2_geo_resonator.gds")
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out_gds_test/")
os.makedirs(OUTDIR, exist_ok=True)

from SQDMetal.PALACE.Frequency_Driven_Simulation import PALACE_Driven_Simulation

# Combined bbox of elements 0,4,5: x[-895.0,107.5] y[-1837.5,2047.5].
# segment_rectangle needs margin beyond the metal bbox so no substrate edge
# coincides exactly with a metal edge (GeomGDS gotcha -- see
# packages/gds-stackup-pipeline README).
MARGIN_UM = 100.0
seg_rect = (-895.0 - MARGIN_UM, -1837.5 - MARGIN_UM, 107.5 + MARGIN_UM, 2047.5 + MARGIN_UM)
seg_rect_m = tuple(v * 1e-6 for v in seg_rect)

driven = PALACE_Driven_Simulation(
    name="gds_phase2_test", sim_parent_directory=OUTDIR,
    mode="simPC", meshing="GMSH", create_files=True,
    gds_design=GDS_PATH,
    dielectric_material="silicon",
    gds_chip_thickness=500e-6,
    segment_rectangle=seg_rect_m,
    user_options={
        "dielectric_material": "silicon", "solns_to_save": 0,
        "solver_order": 2, "solver_tol": 1e-6, "solver_maxits": 500,
        "solver_initial_guess": True, "fillet_resolution": 12,
        "palace_dir": "/dev/null", "num_cpus": 1,
    })
driven.add_metallic(1)
# NOTE: no add_ground_plane() -- GeomGDS ignores it (prints its own warning
# if called), and per this .geo file's own structure there's no separate
# ground-fill element among its 22 pieces anyway (see memory
# sqdmetal-gds-design-ingestion / this session's analysis) -- elements 4/5
# are floating conductors, consistent with the real reference.
driven.set_farfield("absorbing")

# Port at the feedline's simple open tip (top, y=2047.5 -- a clean flat
# two-strip end, unlike the taper/finger cluster at the bottom which this
# crop doesn't attempt to reproduce). Point inset a few um down from the
# true tip, per the "not exactly on the boundary" gotcha.
try:
    driven.create_port_CPW_via_edge_point((0e-6, 2040e-6), 20e-6, impedance_R=50.0)
    driven.set_port_excitation(1)
    port_ok = True
    port_err = None
except Exception as e:
    port_ok = False
    port_err = repr(e)
    print(f"create_port_CPW_via_edge_point FAILED: {port_err}")

MIN_SIZE_UM = 12.0
MAX_SIZE_UM = 190.0
if port_ok:
    driven.fine_mesh_features(MIN_SIZE_UM * 1e-6)
else:
    # Fallback fine-mesh region if the port call itself failed -- still
    # want to see whether meshing alone succeeds on this geometry.
    driven.fine_mesh_in_rectangle(*(v * 1e-6 for v in seg_rect))

driven.set_zBoundary_as_proportion(1.0, 0.3)
driven.set_xBoundary_as_proportion(0.1)
driven.set_yBoundary_as_proportion(0.1)
driven.set_freq_values(freq_start=7.0e9, freq_end=7.3e9, freq_step=0.1e9)

print(f"port_ok={port_ok}")
tp = None
try:
    driven.prepare_simulation()
    print("prepare_simulation() OK")
    mesh_ok = True
    mesh_err = None
except Exception as e:
    import traceback
    mesh_ok = False
    mesh_err = traceback.format_exc()
    print(f"prepare_simulation() FAILED:\n{mesh_err}")

print(f"RESULT: port_ok={port_ok} mesh_ok={mesh_ok}")
