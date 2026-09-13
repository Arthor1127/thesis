"""HFSS *eigenmode* Qc of the hanger readout resonator, with a lumped-RLC
R = 50 Ohm boundary standing in for each 50 Ohm CPW feedline port -- the
follow-up asked for in ../EIGENMODE_QC_REQUEST.md.

Why: the driven S21/S11 sweep on this geometry never resolves the resonance
(|S11| flat at 0.054 +/- 0.002 -- the coupling through the CoupledLineTee's
~15 um ground strip is too weak for a driven adaptive sweep), and the
existing no-port eigenmode run (../ansys_results/testC_eigen.py) has only the
silicon tan_delta = 1.2e-5 as a loss channel, so its Q = 92 231 is a
dielectric-limited Qi, not a coupling Qc.  Adding a resistive boundary at
each feedline port gives the mode a real external loss channel; the complex
eigenfrequency f = fr + i*fi then yields  Qc = fr / (2 fi)  (with the
dielectric tan_delta also present -- separated out numerically below).

Method (matches ../ansys_results/run_ansys_tests.build_from_gds geometry):
  * zero-thickness PEC metal from the GDS, grounded PEC backplane, compact
    box; the hanging resonator meander is in the ground net (it is
    galvanically shorted to ground through res_short), the feedline is the
    signal conductor.
  * at each domain x-face (x = +/- X, X = outer launchpad edge -- the same
    location the driven wave ports use) a lumped-RLC sheet across EACH CPW
    gap, signal -> ground, R = 100 Ohm per gap so the two in parallel
    present 50 Ohm between the feedline conductor and ground (a matched
    ~50 Ohm termination; the field impedance at that face is ~46 Ohm from
    the driven runs).  RES_R env overrides the 50 Ohm effective value.

  conda activate qiskit_qc
  python eigen_qc.py                      # all 3 geometries
  python eigen_qc.py cs22p5_cl250         # one
  RES_NMODES=3 DF=0.05 python eigen_qc.py cs14_cl300
"""
import json
import os
import sys
import time
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
SHARED = os.path.join(ROOT, "ansys_results")
sys.path.insert(0, SHARED)
sys.path.insert(0, HERE)

# Lossless dielectric + PEC metal by default: the ONLY loss channel is then the
# R = 50 port boundary, so the eigenmode's  Q = Re(f)/(2 Im(f))  IS Qc directly
# -- no fragile "subtract the dielectric Qi" step (which is hopeless when the
# coupling is weak and Q_total ~ Qi).  Set SI_TANDELTA=1.2e-5 to reproduce the
# PALACE / no-port loss budget instead.
os.environ.setdefault("SI_TANDELTA", "0.0")

import run_ansys_tests as R          # noqa: E402  shared HFSS helpers
from gds_geometry import read_metal_polygons, bbox   # noqa: E402

GEO = os.path.join(ROOT, "geometry")
CPW_W_UM = 20.0
CPW_G_UM = 11.5

# per-geometry PALACE references from ../EIGENMODE_QC_REQUEST.md
GEOMS = {
    "cs14_cl300": dict(
        gds="full_design_cs14_cl300.gds",
        palace={"eigenmode_Qc_2port_R50": 13562, "eigenmode_Qc_1port_open2nd": 26966,
                "driven_S21_fit_Qc": 15684},
        note="Part 1. PALACE numbers not fully trusted; independent check.",
    ),
    "cs22p5_cl250": dict(
        gds="full_design_cs22p5_cl250.gds",
        palace={"eigenmode_Qc": 95200},
        note="Part 2 PRIMARY. PALACE Qc just past its confirmed ~82-85k ceiling.",
    ),
    "cs25_cl100": dict(
        gds="full_design_cs25_cl100.gds",
        palace={"eigenmode_Qc": 121515},
        note="Part 2 secondary. PALACE Qc well into its saturated zone.",
    ),
}


def _measure(gds_path):
    p = read_metal_polygons(gds_path)
    return {k: bbox(p[k]) for k in ("ground", "pads", "trace")}, p


def _gap_edges_at_face(polys, xf, y_probe_max):
    """(sig_y, gnd_y): signal-conductor outer |y| and ground *inner* |y| at the
    port x-face.  The ground inner edge is the smallest ground vertex |y| that
    lies OUTSIDE the signal conductor -- at the launchpad that is ~275 um
    (pad_gap = 150), not the 11.5 um feedline gap.  A thin x-slice of each
    ground polygon is used so a far-away slot vertex can't be mistaken for the
    local edge."""
    sig_y = CPW_W_UM / 2.0
    for poly in polys["trace"] + polys["pads"]:
        if poly[:, 0].min() - 1 <= xf <= poly[:, 0].max() + 1:
            sig_y = max(sig_y, float(np.max(np.abs(poly[:, 1]))))
    cand = []
    for poly in polys["ground"]:
        if not (poly[:, 0].min() - 1 <= xf <= poly[:, 0].max() + 1):
            continue
        near_x = np.abs(poly[:, 0] - xf) < 80.0
        yy = np.abs(poly[near_x, 1]) if near_x.any() else np.abs(poly[:, 1])
        out = yy[(yy > sig_y + 1.0) & (yy < y_probe_max)]
        if out.size:
            cand.append(float(out.min()))
    gnd_y = min(cand) if cand else sig_y + CPW_G_UM
    return sig_y, gnd_y


def add_feedline_rlc(h, polys, X, Y, r_eff_ohm, strip_um=12.0, over_um=3.0):
    """Two lumped-RLC sheets per port face (one per CPW gap), signal->ground,
    each 2*r_eff so the pair present r_eff between conductor and ground.  Each
    sheet is grown `over_um` into the signal conductor and into the ground so
    the lumped element is galvanically connected to both (an edge exactly
    coincident with a PEC sheet edge is not reliably picked up)."""
    m = h.modeler
    r_gap = 2.0 * r_eff_ohm
    made = []
    for i, xf in enumerate((-X, X), start=1):
        xin = xf - strip_um if xf > 0 else xf + strip_um
        x0, x1 = min(xf, xin), max(xf, xin)
        sig_y, gnd_y = _gap_edges_at_face(polys, xf, 0.9 * Y)
        for side in ("t", "b"):
            nm = f"rlc_p{i}{side}"
            ya, yb = sig_y - over_um, gnd_y + over_um       # overlap both metals
            if side == "b":
                ya, yb = -yb, -ya
            sheet = R._rect(m, x0, x1, ya, yb, 0.0, nm + "_sheet")
            s = 1.0 if side == "t" else -1.0
            xm = 0.5 * (x0 + x1)
            h.assign_lumped_rlc_to_sheet(
                sheet, start_direction=[[xm, s * sig_y, 0.0], [xm, s * gnd_y, 0.0]],
                rlc_type="Parallel", resistance=float(r_gap), name=nm)
            made.append(nm)
    return made


def run_one(key, graphical):
    g = GEOMS[key]
    gds = os.path.join(GEO, g["gds"])
    tag = f"eigen_qc_{key}"
    r_eff = float(os.environ.get("RES_R", "50"))
    nmodes = int(os.environ.get("RES_NMODES", "2"))
    df = float(os.environ.get("DF", "0.1"))            # MaxDeltaFreq, percent
    maxp = int(os.environ.get("MAXPASS", "22"))
    minp = int(os.environ.get("MINPASS", "12"))
    fmin = float(os.environ.get("FMIN", "6"))

    b, polys = _measure(gds)
    # PORT_X: clip the domain here instead of at the outer launchpad edge.
    # Setting it to ~965 drops the whole launchpad (taper + pad), cutting the
    # domain volume ~2x -- which is what makes room in 16 GB RAM for the
    # coupler mesh seed -- and puts the R=50 ports on the clean 11.5 um
    # feedline CPW gap instead of the 150 um launchpad-pad gap.
    X = float(os.environ.get("PORT_X",
                             str(max(abs(b["pads"][0]), abs(b["pads"][2])))))
    ty = max(abs(b["trace"][1]), abs(b["trace"][3]))
    yh = float(os.environ.get("YHALF", str(round(ty + 200, -1))))
    print(f"\n{'='*64}\n  {key}   R_eff={r_eff} Ohm  X={X:.0f}um  y_half={yh:.0f}um "
          f"(resonator |y|={ty:.0f}um)\n  PALACE: {g['palace']}\n{'='*64}", flush=True)

    h_vac = float(os.environ.get("HVAC", "300"))
    if os.environ.get("AEDT_GRPC") == "1":
        # COM transport can wedge after repeated OOM crashes ("InsertDesign:
        # The remote procedure call failed"); gRPC is the fallback.
        R.kill_stale_aedt()
        from ansys.aedt.core import settings as _st
        _st.use_grpc_api = True
        from ansys.aedt.core import Hfss as _Hfss
        h = _Hfss(version=R.AEDT_VERSION, non_graphical=not graphical,
                  new_desktop=True, solution_type="Eigenmode",
                  project=f"ansys_{tag}")
    else:
        h = R.open_hfss(f"ansys_{tag}", graphical, solution_type="Eigenmode")
    try:
        R.build_from_gds(h, gds, port_x_um=X, y_half_um=yh, h_vac_um=h_vac,
                         signal_only_near_axis=True, add_ports=False)
        rlc = add_feedline_rlc(h, polys, X, yh, r_eff)
        print(f"   lumped-RLC boundaries: {rlc}", flush=True)

        # Force fine tets in the CoupledLineTee gap.  The resonator<->feedline
        # coupling is evanescent across a ~15 um strip; the eigenmode adaptive
        # refines where the FIELD is strong (inside the resonator) and leaves
        # that gap coarse, so the resonator mode comes out effectively
        # uncoupled (Qc ~ 1e9 = eigensolver noise floor).  A floating z=0 sheet
        # over the coupler region (no boundary assigned -> electrically
        # transparent) with a length-mesh op constrains the local mesh without
        # adding a solid that would intersect the vacuum/substrate boxes.
        if os.environ.get("SEED_COUPLER", "1") == "1":
            # tight window over the CoupledLineTee overlap only (from the GDS:
            # feedline bottom edge y=-10, resonator coupled-arm top edge
            # y~-47..-67, x roughly [-150,130] -- scales a little with
            # coupling_length but this window covers all 3 points)
            sx0 = float(os.environ.get("SEED_X0", "-200"))
            sx1 = float(os.environ.get("SEED_X1", "170"))
            sy0 = float(os.environ.get("SEED_Y0", "-85"))
            sy1 = float(os.environ.get("SEED_Y1", "18"))
            ml = float(os.environ.get("SEED_MAXLEN", "5.0"))
            mm = h.modeler
            zone = R._rect(mm, sx0, sx1, sy0, sy1, 0.0, "coupler_meshzone")
            try:
                h.mesh.assign_length_mesh([zone], inside_selection=False,
                                          maximum_length=ml,
                                          maximum_elements=120000,
                                          name="seed_coupler")
                print(f"   coupler mesh zone x[{sx0},{sx1}] y[{sy0},{sy1}] "
                      f"maxlen {ml} um", flush=True)
            except Exception as e:  # noqa: BLE001
                print("   coupler seed warn:", repr(e), flush=True)
        R.check_validation(h)
        h.save_project()

        setup = h.create_setup("Eig")
        setup.props["MinimumFrequency"] = f"{fmin}GHz"
        setup.props["NumModes"] = nmodes
        setup.props["MaximumPasses"] = maxp
        setup.props["MinimumPasses"] = minp
        setup.props["MaxDeltaFreq"] = df
        setup.props["MinimumConvergedPasses"] = 2
        setup.props["ConvergeOnRealFreq"] = True
        setup.update()
        print(f"   eigenmode: {nmodes} modes from {fmin} GHz, "
              f"passes {minp}-{maxp}, MaxDeltaFreq {df}%", flush=True)
        t0 = time.time()
        h.analyze_setup("Eig", cores=int(os.environ.get("HFSS_CORES", "6")))
        solve_s = time.time() - t0
        try:                       # persist the solved solution before readback
            h.save_project()       # (readback field-recovery can OOM on this box)
        except Exception as e:  # noqa: BLE001
            print("   post-solve save warn:", repr(e), flush=True)

        tand = float(os.environ["SI_TANDELTA"])
        want_passhist = os.environ.get("READ_PASSHIST", "0") == "1"

        def _modes_for(ssn):
            mm = []
            for i in range(1, nmodes + 1):
                rec = {}
                for q, k in ((f"re(Mode({i}))", "f_hz"),
                             (f"im(Mode({i}))", "f_im_hz")):
                    try:
                        sd = h.post.get_solution_data(expressions=[q],
                                                      setup_sweep_name=ssn)
                        if sd and q in sd.expressions:
                            _, v = sd.get_expression_data(q, "real")
                            rec[k] = [float(x) for x in v]
                    except Exception as e:  # noqa: BLE001
                        rec[k + "_err"] = repr(e)
                mm.append(rec)
            return mm

        # final (converged) modes
        modes = []
        for rec in _modes_for("Eig : LastAdaptive"):
            r = {}
            fh = rec.get("f_hz", [None])[-1] if rec.get("f_hz") else None
            fi = rec.get("f_im_hz", [None])[-1] if rec.get("f_im_hz") else None
            if fh and fi:
                r = {"f_hz": fh, "f_im_hz": fi, "f_ghz": fh / 1e9,
                     "Q": abs(fh / (2 * fi))}
            modes.append(r)

        # per-adaptive-pass history of every mode -> lets us see whether Q has
        # actually converged or is drifting on an eigensolver noise floor.
        # OFF by default: the AdaptivePass query makes HFSS recover fields for
        # every pass and can OOM on a 16 GB box.  READ_PASSHIST=1 to enable,
        # or recover it later from the saved project with read_solved.py.
        pass_hist = []
        if want_passhist:
            try:
                ph = _modes_for("Eig : AdaptivePass")
                npass = max((len(m.get("f_hz", [])) for m in ph), default=0)
                for p in range(npass):
                    row = {"pass": p + 1}
                    for i, m in enumerate(ph, start=1):
                        fh = m.get("f_hz", [])
                        fi = m.get("f_im_hz", [])
                        if p < len(fh) and p < len(fi) and fi[p]:
                            row[f"mode{i}_fghz"] = round(fh[p] / 1e9, 6)
                            row[f"mode{i}_Q"] = round(abs(fh[p] / (2 * fi[p])), 1)
                    pass_hist.append(row)
            except Exception as e:  # noqa: BLE001
                pass_hist = [{"error": repr(e)}]

        # resonator mode: nearest to ~8.3 GHz, but only among modes in a
        # plausible window (7-9.5 GHz) so a spurious low box mode isn't picked
        cand = [m for m in modes if m.get("f_ghz") and 7.0 < m["f_ghz"] < 9.6]
        if not cand:
            cand = [m for m in modes if m.get("f_ghz")]
        res_mode = min(cand, key=lambda m: abs(m["f_ghz"] - 8.3)) if cand else {}

        # Qc.  Lossless (tand==0): Q == Qc directly.  Lossy: subtract the
        # dielectric Qi (only meaningful when Q_total is well below Qi).
        QI_DIEL = 92231.0
        q_tot = res_mode.get("Q")
        if tand == 0.0:
            qc = q_tot
            qc_note = ("dielectric lossless + PEC metal, so the eigenmode Q is "
                       "Qc directly (fr/(2 fi), fi from the R=50 ports only). "
                       "Check pass_history: if mode Q keeps climbing / jumping "
                       "pass-to-pass it has hit the eigensolver noise floor and "
                       "Qc is only a lower bound.")
        else:
            inv = (1.0 / q_tot - 1.0 / QI_DIEL) if q_tot else None
            qc = (1.0 / inv) if (inv and inv > 0) else None
            qc_note = (f"tan_delta={tand}; Qc from 1/Qc = 1/Q_total - 1/Qi, "
                       f"Qi={QI_DIEL}. Unreliable when Q_total ~ Qi.")

        # noise-floor flag: Q of the resonator mode over the last 4 passes
        res_i = modes.index(res_mode) + 1 if res_mode in modes else 1
        qseq = [r.get(f"mode{res_i}_Q") for r in pass_hist[-4:]
                if r.get(f"mode{res_i}_Q")]
        q_converged = None
        if len(qseq) >= 3:
            spread = (max(qseq) - min(qseq)) / max(np.mean(qseq), 1.0)
            q_converged = bool(spread < 0.15)

        out = {
            "geometry": key, "tag": tag, "r_eff_ohm": r_eff,
            "r_per_gap_ohm": 2.0 * r_eff, "n_gaps_per_port": 2, "n_ports": 2,
            "si_tandelta": tand,
            "palace_reference": g["palace"], "geometry_note": g["note"],
            "resonator_mode": {
                "fr_hz": res_mode.get("f_hz"), "fr_ghz": res_mode.get("f_ghz"),
                "fi_hz": res_mode.get("f_im_hz"),
                "Q_total_hfss": q_tot,
                "Qc": qc,
                "Qc_is_lower_bound_only": (q_converged is False),
                "Qc_note": qc_note,
            },
            "pass_history": pass_hist,
            "all_modes": modes, "solve_s": round(solve_s, 1),
        }
        with open(os.path.join(HERE, f"{tag}_result.json"), "w") as fp:
            json.dump(out, fp, indent=2)
        print(json.dumps({k: v for k, v in out.items()
                          if k not in ("all_modes", "pass_history")}, indent=2),
              flush=True)
        for r in pass_hist:
            print("   pass", r, flush=True)
        return out
    finally:
        R.close_hfss(h)


def main():
    sel = [a for a in sys.argv[1:] if not a.startswith("-")] or list(GEOMS)
    graphical = "--graphical" in sys.argv
    allres = {}
    for key in sel:
        try:
            allres[key] = run_one(key, graphical)
        except Exception as e:  # noqa: BLE001
            allres[key] = {"error": repr(e), "traceback": traceback.format_exc()}
            print(allres[key]["traceback"], flush=True)
        with open(os.path.join(HERE, "_eigen_qc_all.json"), "w") as fp:
            json.dump(allres, fp, indent=2, default=str)


if __name__ == "__main__":
    main()
