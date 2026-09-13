"""Eigenmode solve of the Test C hanger-resonator design to pin the resonator
mode frequency f_r (and its unloaded Q from the dielectric loss).  This breaks
the chicken-and-egg of "need f_r to place a narrow driven sweep, need a driven
sweep to find f_r".

No ports: the domain-box PEC walls terminate the feedline (its own modes land
near ~17 GHz, well clear of the ~8 GHz resonator mode).  Silicon tan_delta =
1.2e-5 (SI_TANDELTA env) so the modes carry a finite Q.

  python testC_eigen.py            # 4 modes from 6 GHz
  NMODES=6 FMIN=5 python testC_eigen.py
"""
import json
import os
import sys

from ansys.aedt.core import settings
settings.use_grpc_api = False
from ansys.aedt.core import Hfss

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_ansys_tests as R  # noqa: E402

os.environ.setdefault("SI_TANDELTA", "1.2e-5")
NMODES = int(os.environ.get("NMODES", "4"))
FMIN = float(os.environ.get("FMIN", "6"))
GDS = os.path.join(R.GEO, "full_design_cs14_cl300.gds")

b = R._measure(GDS)
X = max(abs(b["pads"][0]), abs(b["pads"][2]))
ty = max(abs(b["trace"][1]), abs(b["trace"][3]))
yh = round(ty + 250, -1)

h = R.open_hfss("ansys_testC_eig", graphical=True, solution_type="Eigenmode")
try:
    R.build_from_gds(h, GDS, port_x_um=X, y_half_um=yh, h_vac_um=400.0,
                     signal_only_near_axis=True, add_ports=False)
    h.save_project()
    setup = h.create_setup("Eig")
    setup.props["MinimumFrequency"] = f"{FMIN}GHz"
    setup.props["NumModes"] = NMODES
    setup.props["MaximumPasses"] = 18
    setup.props["MaxDeltaFreq"] = 2.0        # percent
    setup.props["ConvergeOnRealFreq"] = True
    setup.update()
    print(f"eigenmode: {NMODES} modes from {FMIN} GHz, y_half={yh}", flush=True)
    h.analyze_setup("Eig", cores=int(os.environ.get("HFSS_CORES", "6")))

    d = h.post.get_solution_data(expressions=["Freq"], setup_sweep_name="Eig : LastAdaptive",
                                 report_category="Eigenmode")
    out = {"modes": []}
    try:
        modes_f = d.primary_sweep_values if d else []
    except Exception:
        modes_f = []
    # more robust: read the eigenmode report table
    try:
        sol = h.post.get_solution_data(expressions=[f"Mode({i})" for i in range(1, NMODES + 1)])
    except Exception:
        sol = None
    for i in range(1, NMODES + 1):
        rec = {}
        for q, key in (("re(Mode(%d))" % i, "f_hz"), ("im(Mode(%d))" % i, "f_imag_hz")):
            try:
                sd = h.post.get_solution_data(expressions=[q],
                                              setup_sweep_name="Eig : LastAdaptive")
                if sd and q in sd.expressions:
                    _, v = sd.get_expression_data(q, "real")
                    rec[key] = float(list(v)[0])
            except Exception as e:
                rec[key + "_err"] = repr(e)
        if rec.get("f_hz") and rec.get("f_imag_hz"):
            rec["Q"] = abs(rec["f_hz"] / (2 * rec["f_imag_hz"]))
        out["modes"].append(rec)
    out["note"] = ("Identify the resonator mode as the one near ~8 GHz whose "
                   "field is concentrated on the hanging meander (feedline "
                   "modes ~17 GHz, box modes >12 GHz).")
    out["palace_reference"] = {"fr_hz": 8.101e9, "eigenmode_fr_range_ghz": [8.25, 8.28]}
    json.dump(out, open(os.path.join(HERE, "testC_eigenmodes.json"), "w"), indent=2)
    print(json.dumps(out, indent=2), flush=True)
finally:
    R.close_hfss(h)
