"""Diagnostic: does a z=0 lumped-RLC sheet across each CPW gap actually
terminate the line?  Wave port at x=-X, my lumped-RLC pair at x=+X, on the
uniform-ish with_launchpads geometry.  If the RLC end is a matched ~50 ohm
load, |S11| should be small (<= -10 dB); if the RLC does nothing, |S11| ~ 0 dB
and the line rings.
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "ansys_results"))
sys.path.insert(0, HERE)
import run_ansys_tests as R                       # noqa: E402
from gds_geometry import read_metal_polygons, bbox  # noqa: E402
import eigen_qc as EQ                              # noqa: E402

GDS = os.path.join(HERE, "..", "geometry", "with_launchpads.gds")
RE = float(os.environ.get("RES_R", "50"))

polys = read_metal_polygons(GDS)
b = {k: bbox(polys[k]) for k in ("ground", "pads", "trace")}
X = max(abs(b["pads"][0]), abs(b["pads"][2]))
Y = 550.0

h = R.open_hfss("ansysport_validation", graphical=False, solution_type="Modal")
try:
    m = h.modeler
    R.build_from_gds(h, GDS, port_x_um=X, y_half_um=Y, h_vac_um=300.0,
                     add_ports=False)

    # wave port at x = -X (full face)
    ps = R._yzrect(m, -X, -Y, Y, -R.SUB_THICK_UM, 300.0, "wpsheet")
    sig_y, gnd_y = EQ._gap_edges_at_face(polys, -X, 0.9 * Y)
    h.wave_port(assignment=ps, integration_line=[[-X, 0.0, 0.0], [-X, gnd_y - 0.5, 0.0]],
                modes=1, renormalize=False, name="p1",
                characteristic_impedance="Zpi")

    # my lumped-RLC pair at x = +X
    rlc = EQ.add_feedline_rlc(h, polys, X, Y, RE)
    print("RLC:", rlc, flush=True)
    R.check_validation(h)
    h.save_project()

    setup = h.create_setup("Setup1")
    setup.props["Frequency"] = "8GHz"
    setup.props["MaximumPasses"] = 12
    setup.props["MaxDeltaS"] = 0.02
    setup.update()
    h.create_linear_count_sweep(setup="Setup1", unit="GHz", start_frequency=6.0,
                                stop_frequency=10.0, num_of_freq_points=201,
                                name="Sweep1", sweep_type="Interpolating",
                                save_fields=False)
    R.setup_and_sweep.last_sweep = "Sweep1"
    h.analyze_setup("Setup1", cores=4)

    import skrf
    tp = os.path.join(HERE, "port_validation.s1p")
    h.export_touchstone(setup="Setup1", sweep="Sweep1", output_file=tp)
    n = skrf.Network(tp)
    s11 = n.s[:, 0, 0]
    db = 20 * np.log10(np.abs(s11))
    print(f"\nRES_R={RE}  |S11| over 6-10 GHz:  mean {db.mean():.1f} dB   "
          f"min {db.min():.1f} dB   max {db.max():.1f} dB", flush=True)
    print("VERDICT:", "RLC ABSORBS (ports work)" if db.mean() < -6
          else "RLC does little (ports NOT terminating)", flush=True)
finally:
    R.close_hfss(h)
