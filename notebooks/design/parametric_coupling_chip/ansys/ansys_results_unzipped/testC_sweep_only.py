"""Reuse Test C's already-solved adaptive mesh (ansys_testC.aedt / design
HFSS_NJY) and run just a sweep on it -- for scanning the resonance without
paying the ~15 min adaptive pass again.

  python testC_sweep_only.py 6.5 8.7 Fast 11001
  python testC_sweep_only.py 7.95 8.15 Discrete 2001
"""
import os
import sys

from ansys.aedt.core import settings
settings.use_grpc_api = False
from ansys.aedt.core import Hfss

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_ansys_tests as R  # noqa: E402

f0 = float(sys.argv[1]) if len(sys.argv) > 1 else 6.5
f1 = float(sys.argv[2]) if len(sys.argv) > 2 else 8.7
kind = sys.argv[3] if len(sys.argv) > 3 else "Fast"
npts = int(sys.argv[4]) if len(sys.argv) > 4 else 11001
tag = os.environ.get("TAG", "testC_full_design") + os.environ.get("SUFFIX", "")

R.kill_stale_aedt()
h = Hfss(version="2024.2", non_graphical=False, new_desktop=True,
         project=r"C:\Users\artho\Documents\Ansoft\ansys_testC.aedt",
         solution_type="Modal")
try:
    sn = h.setup_names[0]
    print("setup:", sn, "existing sweeps:", h.get_sweeps(sn), flush=True)
    for sw in list(h.get_sweeps(sn)):
        try:
            h.get_setup(sn).delete_sweep(sw)
        except Exception as e:
            print("  could not delete", sw, e)
    h.create_linear_count_sweep(setup=sn, unit="GHz", start_frequency=f0,
                                stop_frequency=f1, num_of_freq_points=npts,
                                name="Sweep1", sweep_type=kind, save_fields=False)
    R.setup_and_sweep.last_sweep = "Sweep1"
    print(f"sweep {kind} {f0}-{f1} GHz x{npts} -- solving", flush=True)
    h.analyze_setup(sn, cores=int(os.environ.get("HFSS_CORES", "6")))
    res, data = R.extract(h, tag, HERE, two_port=True)
    R._fit_notch(data, HERE, out_name=tag + "_fit.json")
    import json
    print(json.dumps({k: res[k] for k in res
                      if k not in ("convergence", "analytic_targets")},
                     indent=2, default=str), flush=True)
finally:
    R.close_hfss(h)
