"""HFSS *eigenmode* Qc of the bare shorted quarter-wave CPW resonator, with a
lumped series resistance R across the open-end gap.  Lossless dielectric + PEC
metal, so the only loss channel is R and the eigenmode's complex frequency
gives  Qc = Re(f) / (2 Im(f))  directly.

This is the same method class PALACE used (complex-eigenfrequency), run in a
different tool -- HFSS -- as the cross-check.  It is the fallback for when the
driven S11 sweep cannot resolve the (near-lossless, weakly-loaded) resonance
(see run_bareres_s11.py / SUMMARY.md).

  python bareres_eigen.py            # all three R
  python bareres_eigen.py 1e6
"""
import json
import os
import sys
import time
import traceback

import numpy as np

os.environ["SI_TANDELTA"] = "0.0"     # lossless -> eigenmode Q == Qc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_bareres_s11 as B  # noqa: E402
import run_ansys_tests as R  # noqa: E402  (shared HFSS helpers)

R_VALUES = B.R_VALUES
ANALYTIC_QC = B.ANALYTIC_QC
PALACE_QC = B.PALACE_QC


def run_one(Rohm, outdir, graphical):
    tag = f"bareres_eig_R{Rohm:.0e}".replace("+", "")
    print(f"\n{'='*60}\n  EIGENMODE  R = {Rohm:.3g} Ohm   analytic Qc = "
          f"{ANALYTIC_QC.get(Rohm)}   PALACE {PALACE_QC.get(Rohm)}\n{'='*60}",
          flush=True)
    h = R.open_hfss(f"ansys_{tag}", graphical, solution_type="Eigenmode")
    try:
        geo = B.build(h)                       # same geometry builder
        # shunt resistance R across ONE open-end CPW gap (trace -> ground) --
        # exactly the loading the analytic Qc = (pi/4)(R/Z0) assumes
        try:
            h.assign_lumped_rlc_to_sheet(
                geo["port_sheet"], start_direction=B._iline(geo),
                rlc_type="Parallel", resistance=float(Rohm), name="lumpedR")
        except Exception as e:  # noqa: BLE001
            print("   rlc assign err:", repr(e))
            raise
        h.save_project()

        setup = h.create_setup("Eig")
        setup.props["MinimumFrequency"] = "6GHz"
        setup.props["NumModes"] = 1
        setup.props["MaximumPasses"] = int(os.environ.get("BR_MAXPASS", "26"))
        setup.props["MinimumPasses"] = int(os.environ.get("BR_MINPASS", "16"))
        setup.props["MaxDeltaFreq"] = float(os.environ.get("BR_DF", "0.02"))  # %
        setup.props["MinimumConvergedPasses"] = 2
        setup.props["ConvergeOnRealFreq"] = True
        setup.update()
        t0 = time.time()
        h.analyze_setup("Eig", cores=int(os.environ.get("HFSS_CORES", "6")))
        solve_s = time.time() - t0

        modes = []
        for i in (1,):
            rec = {}
            for q, k in ((f"re(Mode({i}))", "f_hz"), (f"im(Mode({i}))", "f_im_hz")):
                try:
                    sd = h.post.get_solution_data(expressions=[q],
                                                  setup_sweep_name="Eig : LastAdaptive")
                    if sd and q in sd.expressions:
                        _, v = sd.get_expression_data(q, "real")
                        rec[k] = float(list(v)[0])
                except Exception as e:  # noqa: BLE001
                    rec[k + "_err"] = repr(e)
            if rec.get("f_hz") and rec.get("f_im_hz"):
                rec["Q"] = abs(rec["f_hz"] / (2 * rec["f_im_hz"]))
            modes.append(rec)
        m0 = modes[0]
        res = {
            "R_ohm": Rohm, "tag": tag, "method": "HFSS eigenmode, lossless, lumped series R",
            "analytic_Qc": ANALYTIC_QC.get(Rohm),
            "analytic_formula": "Qc = (pi/4)*(R/49.961)",
            "palace_eigenmode_Qc": PALACE_QC.get(Rohm),
            "fr_ghz": round(m0.get("f_hz", 0) / 1e9, 5) if m0.get("f_hz") else None,
            "Qc_eigenmode": round(m0["Q"], 1) if m0.get("Q") else None,
            "Qc_over_analytic": (round(m0["Q"] / ANALYTIC_QC[Rohm], 4)
                                 if m0.get("Q") and ANALYTIC_QC.get(Rohm) else None),
            "modes": modes, "solve_s": round(solve_s, 1),
        }
        with open(os.path.join(outdir, f"{tag}_result.json"), "w") as fp:
            json.dump(res, fp, indent=2)
        print(json.dumps({k: v for k, v in res.items() if k != "modes"}, indent=2),
              flush=True)
        return res
    finally:
        R.close_hfss(h)


def main():
    rs = [float(x) for x in sys.argv[1:] if not x.startswith("-")] or R_VALUES
    allres = {}
    for Rohm in rs:
        try:
            allres[f"{Rohm:.0e}"] = run_one(Rohm, HERE, True)
        except Exception as e:  # noqa: BLE001
            allres[f"{Rohm:.0e}"] = {"error": repr(e),
                                     "traceback": traceback.format_exc()}
            print(allres[f"{Rohm:.0e}"]["traceback"], flush=True)
        with open(os.path.join(HERE, "_eig_all_results.json"), "w") as fp:
            json.dump(allres, fp, indent=2, default=str)


if __name__ == "__main__":
    main()
