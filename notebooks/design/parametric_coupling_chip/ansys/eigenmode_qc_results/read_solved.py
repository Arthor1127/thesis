"""Recover the eigenmode complex frequencies from an already-SOLVED
eigen_qc project, when eigen_qc.py's own readback crashed HFSS (field
recovery OOM on this 16 GB box).  Opens the project read-only, no re-solve.

  conda activate qiskit_qc
  python read_solved.py cs14_cl300 [nmodes]
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "ansys_results"))
sys.path.insert(0, HERE)
import run_ansys_tests as R          # noqa: E402
import eigen_qc as EQ                # noqa: E402

key = sys.argv[1] if len(sys.argv) > 1 else "cs14_cl300"
nmodes = int(sys.argv[2]) if len(sys.argv) > 2 else 3
tag = f"eigen_qc_{key}"
proj = os.path.expanduser(f"~/Documents/Ansoft/ansys_{tag}.aedt")

R.kill_stale_aedt()
from ansys.aedt.core import settings
settings.use_grpc_api = False
from ansys.aedt.core import Hfss

h = Hfss(version="2024.2", non_graphical=True, new_desktop=True,
         project=proj, solution_type="Eigenmode")
try:
    modes = []
    for i in range(1, nmodes + 1):
        rec = {}
        for q, k in ((f"re(Mode({i}))", "f_hz"), (f"im(Mode({i}))", "f_im_hz")):
            try:
                sd = h.post.get_solution_data(
                    expressions=[q], setup_sweep_name="Eig : LastAdaptive")
                if sd and q in sd.expressions:
                    _, v = sd.get_expression_data(q, "real")
                    rec[k] = float(list(v)[-1])
            except Exception as e:  # noqa: BLE001
                rec[k + "_err"] = repr(e)
        if rec.get("f_hz") and rec.get("f_im_hz"):
            rec["f_ghz"] = rec["f_hz"] / 1e9
            rec["Q"] = abs(rec["f_hz"] / (2 * rec["f_im_hz"]))
        modes.append(rec)
    print(json.dumps(modes, indent=2), flush=True)

    cand = [m for m in modes if m.get("f_ghz") and 7.0 < m["f_ghz"] < 9.6] \
        or [m for m in modes if m.get("f_ghz")]
    res = min(cand, key=lambda m: abs(m["f_ghz"] - 8.3)) if cand else {}

    g = EQ.GEOMS[key]
    tand = float(os.environ.get("SI_TANDELTA", "0.0"))
    q_tot = res.get("Q")
    QI = 92231.0
    if tand == 0.0:
        qc = q_tot
    else:
        inv = (1.0 / q_tot - 1.0 / QI) if q_tot else None
        qc = (1.0 / inv) if (inv and inv > 0) else None
    out = {
        "geometry": key, "tag": tag, "si_tandelta": tand,
        "palace_reference": g["palace"], "geometry_note": g["note"],
        "recovered_from": "read_solved.py (eigen_qc readback had crashed)",
        "resonator_mode": {
            "fr_hz": res.get("f_hz"), "fr_ghz": res.get("f_ghz"),
            "fi_hz": res.get("f_im_hz"), "Q_total_hfss": q_tot, "Qc": qc,
        },
        "all_modes": modes,
    }
    with open(os.path.join(HERE, f"{tag}_result.json"), "w") as fp:
        json.dump(out, fp, indent=2)
    print("wrote", f"{tag}_result.json", flush=True)
finally:
    R.close_hfss(h)
