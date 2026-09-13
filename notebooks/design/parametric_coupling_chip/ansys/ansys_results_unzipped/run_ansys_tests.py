"""ANSYS HFSS cross-check of the qiskit_metal / SQDMetal / PALACE CPW impedance
discrepancy -- see ../HANDOFF.md.

Builds a clean 3D HFSS model (zero-thickness PEC metal on a silicon slab with a
grounded backplane) directly from the handoff GDS metal polygons, terminates
each end with an HFSS *wave port* (whose characteristic impedance Zpi comes
straight from the 2D port field solution, an independent methodology from
PALACE's lumped-port S-parameter Zin), runs a driven-modal solve, and writes
the raw results into this directory.

  Test A  line_only              bare CPW, wave ports at the trace ends
  Test B  with_launchpads        + launchpad tapers, ports at the outer pad edge
  Test C  full_design_cs14_cl300 full hanger-resonator design, 2-port S21

Usage (in the qiskit_qc conda env):
  python run_ansys_tests.py A            # one test
  python run_ansys_tests.py A B
  python run_ansys_tests.py all
Options:
  --graphical     show the AEDT window (default: shown; COM needs a session)
  --no-solve      build + validate + save only, skip the solve

Env:  AEDT_VERSION (default 2024.2)   HFSS_CORES (default 4)
"""
import argparse
import json
import os
import sys
import time
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
GEO = os.path.abspath(os.path.join(HERE, "..", "geometry"))
sys.path.insert(0, HERE)
from gds_geometry import read_metal_polygons, bbox  # noqa: E402

AEDT_VERSION = os.environ.get("AEDT_VERSION", "2024.2")
CORES = int(os.environ.get("HFSS_CORES", "4"))
QUICK = False  # set by --quick: tiny mesh + few freq points, for pipeline shakeout
ZDEF = os.environ.get("HFSS_ZDEF", "Zpi")  # wave-port char.-impedance definition

# cross-section / material -- ../our_results/impedance_results_summary.json
ER_SILICON = 11.45
SUB_THICK_UM = 500.0
CPW_W_UM = 20.0
CPW_G_UM = 11.5

# ---------------------------------------------------------------------------
# AEDT connection (COM / PythonNET -- gRPC does not connect on this machine)
# ---------------------------------------------------------------------------

def kill_stale_aedt():
    """pyaedt's COM release on this machine often leaves the ansysedt.exe it
    spawned running; a stale process then breaks the next COM launch
    (IndexError in __init_dotnet). Clear them before each run."""
    os.system('taskkill /F /IM ansysedt.exe >NUL 2>&1')
    os.system('taskkill /F /IM ansyscl.exe >NUL 2>&1')
    time.sleep(1.0)
    import glob
    for lk in glob.glob(os.path.expanduser("~/Documents/Ansoft/ansys_test*.aedt.lock")):
        try:
            os.remove(lk)
        except OSError:
            pass


def open_hfss(project, graphical=True, solution_type="Modal"):
    kill_stale_aedt()
    from ansys.aedt.core import settings
    settings.use_grpc_api = False          # gRPC does not connect on this machine
    from ansys.aedt.core import Hfss
    h = Hfss(version=AEDT_VERSION, non_graphical=not graphical, new_desktop=True,
             solution_type=solution_type, project=project)
    return h


def close_hfss(h):
    """Release the desktop and make sure the ansysedt process is really gone
    (pyaedt's COM release leaves a zombie on this machine)."""
    pid = None
    try:
        pid = h.desktop_class.aedt_process_id
    except Exception:
        pass
    try:
        h.save_project()
    except Exception:
        pass
    try:
        h.release_desktop(close_projects=True, close_desktop=True)
    except Exception:
        pass
    if pid:
        os.system(f"taskkill /F /PID {pid} >NUL 2>&1")


def check_validation(h):
    msgs, ok = h.validate_full_design()
    for m in msgs:
        print("   [validate]", m)
    return ok


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------

def _p3(poly_um, z=0.0):
    return [[float(x), float(y), float(z)] for x, y in poly_um]


def _rect(m, x0, x1, y0, y1, z, name):
    return m.create_polyline(
        points=[[x0, y0, z], [x1, y0, z], [x1, y1, z], [x0, y1, z]],
        name=name, cover_surface=True, close_surface=True).name


def _yzrect(m, xf, y0, y1, z0, z1, name):
    return m.create_polyline(
        points=[[xf, y0, z0], [xf, y1, z0], [xf, y1, z1], [xf, y0, z1]],
        name=name, cover_surface=True, close_surface=True).name


def _ring_sheet(m, ring, name):
    pts = [[float(x), float(y), 0.0] for x, y in list(ring.coords)[:-1]]
    return m.create_polyline(points=pts, name=name, cover_surface=True,
                             close_surface=True).name


def _draw_shapely(m, geom, prefix):
    """Draw a shapely (Multi)Polygon (holes allowed) as z=0 HFSS sheets.
    Returns the list of resulting sheet names."""
    from gds_geometry import as_polygon_list
    out = []
    for j, poly in enumerate(as_polygon_list(geom)):
        ext = _ring_sheet(m, poly.exterior, f"{prefix}_{j}")
        holes = [_ring_sheet(m, r, f"{prefix}_{j}_h{k}")
                 for k, r in enumerate(poly.interiors)]
        if holes:
            m.subtract([ext], holes, keep_originals=False)
        out.append(ext)
    return out


def build_uniform_cpw(h, length_um, w_um, g_um, y_half_um, h_vac_um, name):
    """Test A: rebuild the bare CPW as clean rectangles (uniform cross-section,
    so the GDS polygon detail is not needed -- w/g/length are taken from it)."""
    m = h.modeler
    m.model_units = "um"
    L = length_um
    xlo, xhi = -L / 2.0, L / 2.0
    yt = w_um / 2.0
    yg = w_um / 2.0 + g_um
    Y = y_half_um

    sil = _silicon(h)
    sub = m.create_box([xlo, -Y, -SUB_THICK_UM], [L, 2 * Y, SUB_THICK_UM],
                       name="substrate", material=sil)
    vac = m.create_box([xlo, -Y, 0.0], [L, 2 * Y, h_vac_um], name="vacuum",
                       material="vacuum")

    trace = _rect(m, xlo, xhi, -yt, yt, 0.0, "trace")
    gnd_l = _rect(m, xlo, xhi, -Y, -yg, 0.0, "gnd_l")
    gnd_r = _rect(m, xlo, xhi, yg, Y, 0.0, "gnd_r")

    h.assign_perfecte_to_sheets([trace], name="PerfE_trace")
    h.assign_perfecte_to_sheets([gnd_l, gnd_r], name="PerfE_gnd")
    # grounded backplane (matches PALACE z_prop_bottom=0 and the CpwParams formula)
    h.assign_perfecte_to_sheets([sub.bottom_face_z.id], name="PerfE_backplane")

    _seed_cpw_mesh(h, [trace, gnd_l, gnd_r])

    ports = []
    for i, xf in enumerate((xlo, xhi), start=1):
        ps = _yzrect(m, xf, -Y, Y, -SUB_THICK_UM, h_vac_um, f"portsheet{i}")
        iline = [[xf, 0.0, 0.0], [xf, yg, 0.0]]
        p = h.wave_port(assignment=ps, integration_line=iline, modes=1,
                        renormalize=False, name=f"p{i}",
                        characteristic_impedance=ZDEF)
        ports.append(p)
    return ports


def _seed_cpw_mesh(h, sheets, max_len=5.0):
    """(kept for reference) A blanket length-mesh on the full-length CPW sheets
    exploded to ~500k tets. HFSS's port-adaptive + lambda refinement already
    converges this geometry (10 passes, ~8k tets, dS<0.02) with Zpi flat to the
    4th figure, so seeding is disabled by default. Enable with SEED_MESH=1 for a
    convergence cross-check on a short run."""
    if os.environ.get("SEED_MESH") != "1":
        return
    try:
        h.mesh.assign_length_mesh(sheets, inside_selection=False,
                                  maximum_length=max_len, maximum_elements=200000,
                                  name="seed_cpw")
    except Exception as e:  # noqa: BLE001
        print("   mesh seed warning:", repr(e))


def build_from_gds(h, gds_path, port_x_um, y_half_um, h_vac_um,
                   signal_only_near_axis=False, axis_band_um=40.0, add_ports=True):
    """Tests B / C: use the real GDS metal polygons, clipped in x to the domain
    (which ends at +/- port_x_um). Wave ports on the two x end faces
    (add_ports=False for an eigenmode solve -- the domain walls terminate the
    feedline instead).

    signal_only_near_axis: for the hanger-resonator design (C), only the
    feedline polygons (|y| within axis_band_um of the CPW axis) are the driven
    signal conductor; the hanging resonator meander is galvanically tied to
    ground through its ShortToGround, so its polygons go into the ground net.
    """
    from shapely.geometry import box as _box
    from gds_geometry import heal, as_polygon_list

    m = h.modeler
    m.model_units = "um"
    polys = read_metal_polygons(gds_path)
    X = float(port_x_um)
    Y = float(y_half_um)
    domain = _box(-X, -Y, X, Y)

    sil = _silicon(h)
    sub = m.create_box([-X, -Y, -SUB_THICK_UM], [2 * X, 2 * Y, SUB_THICK_UM],
                       name="substrate", material=sil)
    vac = m.create_box([-X, -Y, 0.0], [2 * X, 2 * Y, h_vac_um], name="vacuum",
                       material="vacuum")

    # classify trace/pad polygons: for the hanger design (C) the feedline sits on
    # the CPW axis while the hanging resonator meander is entirely off to y<0 and
    # is tied to ground through its ShortToGround.
    sig_raw, gnd_raw = list(polys["pads"]), list(polys["ground"])
    for poly in polys["trace"]:
        ymax, ymin = float(poly[:, 1].max()), float(poly[:, 1].min())
        on_axis = (ymin < axis_band_um and ymax > -axis_band_um
                   and abs(0.5 * (ymax + ymin)) < axis_band_um)
        (sig_raw if (not signal_only_near_axis or on_axis) else gnd_raw).append(poly)

    sig_geom = heal(sig_raw).intersection(domain)
    gnd_geom = heal(gnd_raw).intersection(domain)

    sig_names = _draw_shapely(m, sig_geom, "sig")
    gnd_names = _draw_shapely(m, gnd_geom, "gnd")
    sig = m.unite(sig_names) if len(sig_names) > 1 else sig_names[0]
    gnd = m.unite(gnd_names) if len(gnd_names) > 1 else gnd_names[0]
    h.assign_perfecte_to_sheets([sig], name="PerfE_sig")
    h.assign_perfecte_to_sheets([gnd], name="PerfE_gnd")
    h.assign_perfecte_to_sheets([sub.bottom_face_z.id], name="PerfE_backplane")
    _seed_cpw_mesh(h, [sig, gnd], max_len=6.0)

    # optional targeted mesh refinement on the coupler+resonator region -- the
    # weak hanger coupling is invisible to the driven adaptive without it
    seed_box = os.environ.get("TESTC_SEED_BOX")   # "x0,x1,y0,y1,zdepth,maxlen"
    if seed_box:
        x0, x1, ys0, ys1, zd, ml = [float(v) for v in seed_box.split(",")]
        try:
            sb = m.create_box([x0, ys0, -zd], [x1 - x0, ys1 - ys0, 2 * zd],
                              name="mesh_seed")
            sb.model = False
            if m["mesh_seed"].model:            # setter did not take -> force it
                m.oeditor.ChangeProperty(
                    ["NAME:AllTabs", ["NAME:Geometry3DAttributeTab",
                     ["NAME:PropServers", "mesh_seed"],
                     ["NAME:ChangedProps", ["NAME:Model", "Value:=", False]]]])
            h.mesh.assign_length_mesh(["mesh_seed"], inside_selection=True,
                                      maximum_length=ml, maximum_elements=300000,
                                      name="seed_coupler")
            print(f"   mesh seed box {seed_box} (model={m['mesh_seed'].model})",
                  flush=True)
        except Exception as e:  # noqa: BLE001
            print("   seed box warning:", repr(e), flush=True)

    sig_polys = polys["trace"] + polys["pads"]
    if not add_ports:
        return []
    ports = []
    for i, xf in enumerate((-X, X), start=1):
        # integration line: from the signal centre out to the ground inner edge
        # at this port face (so Vpi = int E.dl spans the whole gap).
        shy = CPW_W_UM / 2.0
        for poly in sig_polys:
            if poly[:, 0].min() - 1 <= xf <= poly[:, 0].max() + 1:
                shy = max(shy, float(np.max(np.abs(poly[:, 1]))))
        gy = shy + CPW_G_UM
        for poly in polys["ground"]:
            if poly[:, 0].min() - 1 <= xf <= poly[:, 0].max() + 1:
                yy = np.abs(poly[:, 1])
                near = yy[yy > shy]
                if near.size:
                    gy = min(gy, float(near.min()))
        y_end = min(gy - 0.5, 0.9 * Y)
        ps = _yzrect(m, xf, -Y, Y, -SUB_THICK_UM, h_vac_um, f"portsheet{i}")
        iline = [[xf, 0.0, 0.0], [xf, y_end, 0.0]]
        p = h.wave_port(assignment=ps, integration_line=iline, modes=1,
                        renormalize=False, name=f"p{i}",
                        characteristic_impedance=ZDEF)
        ports.append(p)
    return ports


def _silicon(h):
    # loss-free for Tests A/B (handoff: safe to treat tan_delta as 0); Test C
    # sets SI_TANDELTA=1.2e-5 to match the PALACE sim so the lossless hanger
    # actually produces an S21 notch rather than a pure phase swing.
    tand = float(os.environ.get("SI_TANDELTA", "0.0"))
    name = "silicon_er1145" if tand == 0 else "silicon_er1145_lossy"
    if name not in h.materials.material_keys:
        mat = h.materials.add_material(name)
        mat.permittivity = ER_SILICON
        mat.dielectric_loss_tangent = tand
    return name


# ---------------------------------------------------------------------------
# solve + extract
# ---------------------------------------------------------------------------

def setup_and_sweep(h, f_solve_ghz, f0_ghz, f1_ghz, npts, kind="Interpolating",
                    max_passes=25, delta_s=0.02):
    if QUICK:
        max_passes, delta_s, npts, kind = 3, 0.05, 5, "Discrete"
    if "Setup1" in (h.setup_names or ()):
        setup = h.get_setup("Setup1")
        setup.props["MaximumPasses"] = max_passes  # keep adapted mesh; 1 = reuse
    else:
        setup = h.create_setup("Setup1")
        setup.props["Frequency"] = f"{f_solve_ghz}GHz"
        setup.props["MaximumPasses"] = max_passes
        setup.props["MaxDeltaS"] = delta_s
        setup.props["MinimumConvergedPasses"] = 2
    setup.update()
    # unique sweep name so a stage-2 sweep doesn't collide with stage 1
    swname = "Sweep1"
    n = 1
    while swname in [s for s in (h.get_sweeps("Setup1") or [])]:
        n += 1
        swname = f"Sweep{n}"
    h.create_linear_count_sweep(setup="Setup1", unit="GHz",
                                start_frequency=f0_ghz, stop_frequency=f1_ghz,
                                num_of_freq_points=npts, name=swname,
                                sweep_type=kind, save_fields=False)
    setup_and_sweep.last_sweep = swname
    return setup


def solve(h):
    t0 = time.time()
    h.analyze_setup("Setup1", cores=CORES)
    return time.time() - t0


def _report_cplx(h, expr, ssn):
    """(freq_hz, complex array) for one Modal-Solution-Data report expression,
    or (None, None)."""
    d = h.post.get_solution_data(expressions=[expr], setup_sweep_name=ssn,
                                 domain="Sweep")
    if not d or expr not in d.expressions:
        return None, None
    x, re = d.get_expression_data(expr, formula="real")
    _, im = d.get_expression_data(expr, formula="imag")
    re = np.asarray(re, float)
    im = np.asarray(im, float)
    if not re.size or not np.isfinite(re).any():
        return None, None
    return np.asarray(x, float) * 1e9, re + 1j * im


def _port_zo(h, res, two_port):
    """Field-based port characteristic impedance Zo(p1[/p2]) and the mode's
    Gamma / Epsilon (eps_eff) vs freq, from the sweep solution."""
    ssn = f"Setup1 : {getattr(setup_and_sweep, 'last_sweep', 'Sweep1')}"
    out = {}
    faxis = None
    names = ["p1", "p2"] if two_port else ["p1"]
    for pn in names:
        for expr in (f"Zo({pn})", f"Zo({pn}:1)"):
            f, z = _report_cplx(h, expr, ssn)
            if z is not None:
                out[pn], faxis = z, f
                break
        else:
            out[pn] = None
            res.setdefault("warnings", []).append(f"Zo({pn}) not available")
    for extra in ("Gamma(p1)", "Epsilon(p1)"):
        f, v = _report_cplx(h, extra, ssn)
        if v is not None:
            out[extra] = v
    return out, faxis


def extract(h, name, outdir, two_port):
    import skrf
    res = {"name": name, "two_port": two_port}

    # --- field-based port Zo (= Zpi, from the 2D port field solution) ---
    zo, zf = _port_zo(h, res, two_port)
    z1 = zo.get("p1")
    z2 = zo.get("p2")

    # --- S-parameters via Touchstone.  HFSS exports these NOT renormalized, i.e.
    # referenced to each port's own Zpi(f); skrf assumes 50 from the header, so we
    # set the true per-frequency reference impedance before renormalizing to a
    # flat 50 ohm (matching PALACE's Z0_ref = 50 convention). ---
    ext = "s2p" if two_port else "s1p"
    tpath = os.path.join(outdir, f"{name}.{ext}")
    _sw = getattr(setup_and_sweep, "last_sweep", "Sweep1")
    h.export_touchstone(setup="Setup1", sweep=_sw, output_file=tpath)
    ntwk = skrf.Network(tpath)
    freqs = ntwk.f.astype(float)
    S11_raw = ntwk.s[:, 0, 0]                     # ref = Zpi(f), as solved

    def _deglitch(z):
        """The HFSS 'Port Zo' report can spike to a few / >100 ohm at a single
        frequency right on a resonance -- renormalising through such a spike
        injects a fake S-parameter notch. Clamp isolated outliers to the local
        median."""
        z = np.asarray(z, complex)
        med = np.median(z.real)
        bad = np.abs(z.real - med) > 0.35 * med
        if bad.any():
            z = z.copy()
            z[bad] = med + 0j
        return z

    if z1 is not None and len(z1) == len(freqs):
        z1c = _deglitch(z1)
        z2c = _deglitch(z2) if (two_port and z2 is not None
                                and len(z2) == len(freqs)) else z1c
        ntwk.z0 = np.column_stack([z1c, z2c])[:, :ntwk.nports]
        ntwk.renormalize(50.0)
    S11 = ntwk.s[:, 0, 0]
    S21 = ntwk.s[:, 1, 0] if two_port else np.full_like(S11, np.nan)
    Zin = 50.0 * (1 + S11) / (1 - S11)
    try:
        ntwk.write_touchstone(os.path.join(outdir, f"{name}_50ohm.{ext}"),
                              write_z0=True)
    except Exception:  # noqa: BLE001
        pass

    np.savetxt(os.path.join(outdir, f"{name}_sparams.csv"),
               np.column_stack([freqs,
                                20 * np.log10(np.abs(S11)), np.degrees(np.angle(S11)),
                                20 * np.log10(np.abs(S21)), np.degrees(np.angle(S21)),
                                np.abs(S11_raw)]),
               delimiter=",",
               header=("freq_hz,S11_mag_db_50ohm,S11_phase_deg_50ohm,"
                       "S21_mag_db_50ohm,S21_phase_deg_50ohm,S11_mag_raw_Zoref"),
               comments="")

    if z1 is not None:
        zf_out = zf if (zf is not None and len(zf) == len(freqs)) else freqs
        cols = [zf_out, z1.real, z1.imag]
        hdr = "freq_hz,re_Zpi_port1,im_Zpi_port1"
        if len(zf_out) == len(freqs):
            cols += [Zin.real, Zin.imag]
            hdr += ",re_Zin_fromS11_50,im_Zin_fromS11_50"
        np.savetxt(os.path.join(outdir, f"{name}_Zo.csv"),
                   np.column_stack(cols), delimiter=",", header=hdr, comments="")

    # eps_eff from the port mode (Epsilon(p1)) -- directly comparable to analytic
    eps = zo.get("Epsilon(p1)")
    gam = zo.get("Gamma(p1)")

    # --- convergence (per adaptive pass) ---
    conv = _convergence(h, name, outdir, res)

    # --- summary numbers ---
    if two_port:
        m = (freqs >= 7.9e9) & (freqs <= 8.3e9)
    else:
        m = (freqs >= 6e9) & (freqs <= 10e9)
    if not m.any():
        m = np.ones_like(freqs, dtype=bool)
    res.update({
        "n_freq": int(len(freqs)),
        "f_min_hz": float(freqs.min()), "f_max_hz": float(freqs.max()),
        "Zpi_port1_mean_real": _m(z1.real if z1 is not None else None),
        "Zpi_port1_mean_imag": _m(z1.imag if z1 is not None else None),
        "Zpi_port1_std_real": _s(z1.real if z1 is not None else None),
        "Zpi_port2_mean_real": _m(zo["p2"].real) if zo.get("p2") is not None else None,
        "eps_eff_port1_mean": _m(eps.real) if eps is not None else None,
        "Zin_fromS11_50_mean_real": float(np.mean(Zin.real[m])),
        "Zin_fromS11_50_mean_imag": float(np.mean(Zin.imag[m])),
        "Zin_fromS11_50_abs_mean": float(np.mean(np.abs(Zin[m]))),
        "S11_50_mean_db": float(np.mean(20 * np.log10(np.abs(S11[m])))),
        "abs_S11_50ohm_mean": float(np.mean(np.abs(S11[m]))),
        "S11_raw_mean_abs": float(np.mean(np.abs(S11_raw[m]))),
        "convergence": conv,
        "analytic_targets": {"DesignHelper_Z0": 49.96, "SQDMetal_CpwParams_Z0": 50.29,
                             "PALACE_driven_Zin_abs": 40.48, "PALACE_cap_Z0": 47.55},
    })
    with open(os.path.join(outdir, f"{name}_result.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res, dict(freqs=freqs, S11=S11, S21=S21, Zpi1=z1)


def _convergence(h, name, outdir, res):
    """Adaptive-pass history: tets solved + max |dS| per pass."""
    out = {}
    try:
        cpath = os.path.join(outdir, f"{name}_convergence.csv")
        h.export_convergence("Setup1", output_file=cpath)
        out["file"] = os.path.basename(cpath)
        with open(cpath) as f:
            out["raw"] = f.read().strip().splitlines()[-12:]
    except Exception as e:  # noqa: BLE001
        res.setdefault("warnings", []).append(f"convergence: {e!r}")
    return out


def _m(a):
    return None if a is None else float(np.mean(a))


def _s(a):
    return None if a is None else float(np.std(a))


# ---------------------------------------------------------------------------
# per-test drivers
# ---------------------------------------------------------------------------

def _measure(gds_path):
    p = read_metal_polygons(gds_path)
    return {k: bbox(p[k]) for k in ("ground", "pads", "trace")}


def run_A(outdir, solve_it, graphical):
    gds = os.path.join(GEO, "line_only.gds")
    b = _measure(gds)
    L = b["trace"][2] - b["trace"][0]
    w = b["trace"][3] - b["trace"][1]
    yh = float(os.environ.get("TESTA_YHALF", "200"))
    hv = float(os.environ.get("TESTA_HVAC", str(min(yh, 400.0))))
    tag = os.environ.get("TESTA_TAG", "testA_line_only")
    gnd_full = b["ground"][3] - b["ground"][1]
    print(f"[A] GDS: trace L={L:.1f}um w={w:.1f}um gap {CPW_G_UM}um; "
          f"y_half={yh}um h_vac={hv}um (GDS ground half-width {gnd_full/2:.0f}um)  tag={tag}")
    h = open_hfss("ansys_testA", graphical)
    try:
        build_uniform_cpw(h, L, w, CPW_G_UM, y_half_um=yh, h_vac_um=hv, name="testA")
        check_validation(h)
        h.save_project()
        if not solve_it:
            return {"name": tag, "built_only": True}
        setup_and_sweep(h, 8.0, 6.0, 10.0, 41, kind="Interpolating")
        dt = solve(h)
        res, _ = extract(h, tag, outdir, two_port=True)
        res["solve_s"] = dt
        res["y_half_um"] = yh
        return res
    finally:
        close_hfss(h)


def run_B(outdir, solve_it, graphical):
    gds = os.path.join(GEO, "with_launchpads.gds")
    b = _measure(gds)
    X = max(abs(b["pads"][0]), abs(b["pads"][2]))
    print(f"[B] GDS: pad outer |x|={X:.1f}um")
    h = open_hfss("ansys_testB", graphical)
    try:
        build_from_gds(h, gds, port_x_um=X, y_half_um=550.0, h_vac_um=300.0)
        check_validation(h)
        h.save_project()
        if not solve_it:
            return {"name": "testB", "built_only": True}
        setup_and_sweep(h, 8.0, 6.0, 10.0, 41, kind="Interpolating")
        dt = solve(h)
        res, _ = extract(h, "testB_launchpads", outdir, two_port=True)
        res["solve_s"] = dt
        return res
    finally:
        close_hfss(h)


def run_C(outdir, solve_it, graphical):
    gds = os.path.join(GEO, "full_design_cs14_cl300.gds")
    b = _measure(gds)
    X = max(abs(b["pads"][0]), abs(b["pads"][2]))
    # the hanging resonator meander reaches y ~ -1317um -- the domain MUST contain
    # it or the lambda/4 resonator gets clipped and stops resonating (that is why
    # earlier Test C runs showed no notch). Size y_half to the full trace extent
    # plus margin.
    ty = max(abs(b["trace"][1]), abs(b["trace"][3]))
    yh = float(os.environ.get("TESTC_YHALF", str(round(ty + 250, -1))))
    print(f"[C] GDS: pad |x|={X:.0f}um  resonator reaches |y|={ty:.0f}um -> "
          f"y_half={yh:.0f}um  SI_TANDELTA={os.environ.get('SI_TANDELTA','0')}")
    h = open_hfss("ansys_testC", graphical)
    try:
        build_from_gds(h, gds, port_x_um=X, y_half_um=yh, h_vac_um=400.0,
                       signal_only_near_axis=True)
        check_validation(h)
        h.save_project()
        if not solve_it:
            return {"name": "testC", "built_only": True}
        tag = os.environ.get("TESTC_TAG", "testC_full_design")
        f0 = float(os.environ.get("TESTC_F0", "7.6"))
        f1 = float(os.environ.get("TESTC_F1", "8.7"))
        fs = float(os.environ.get("TESTC_FSOLVE", "8.15"))
        kind = os.environ.get("TESTC_SWEEP", "Fast")   # ALPS/Pade ROM -- like
        npts = int(os.environ.get("TESTC_NPTS", "5501"))  # PALACE's PROM sweep,
        #                                                  catches sharp poles
        setup_and_sweep(h, fs, f0, f1, npts, kind=kind, max_passes=30)
        dt = solve(h)
        res, data = extract(h, tag, outdir, two_port=True)
        res["solve_s"] = dt
        _fit_notch(data, outdir)

        # Stage 2: if a resonance is visible, confirm with a narrow discrete
        # high-resolution sweep centred on it.
        if os.environ.get("TESTC_NO_STAGE2") != "1":
            fr = _locate_feature(data)
            if fr is not None:
                lo, hi = (fr - 8e6) / 1e9, (fr + 8e6) / 1e9
                print(f"[C] stage-2 discrete confirm {lo:.5f}-{hi:.5f} GHz", flush=True)
                setup_and_sweep(h, fr / 1e9, lo, hi, 1601, kind="Discrete",
                                max_passes=1)
                solve(h)
                res2, data2 = extract(h, tag + "_zoom", outdir, two_port=True)
                _fit_notch(data2, outdir, out_name=tag + "_zoom_fit.json")
                res["stage2_window_ghz"] = [lo, hi]
        return res
    finally:
        close_hfss(h)


def _locate_feature(data):
    """Frequency (Hz) of the sharpest S21 amplitude/phase feature, or None if the
    band is flat (span < 0.2 dB and < 20 deg of excess phase)."""
    f = np.asarray(data["freqs"], float)
    s = np.asarray(data["S21"])
    if s.size < 16 or np.any(~np.isfinite(s)):
        return None
    mag_db = 20 * np.log10(np.abs(s))
    ph = np.unwrap(np.angle(s))
    ph_detrend = ph - np.polyval(np.polyfit(f, ph, 1), f)
    if (mag_db.max() - mag_db.min() < 0.2
            and np.degrees(ph_detrend.max() - ph_detrend.min()) < 20):
        return None
    # sharpest point = max local curvature of |S21| dB or detrended phase
    score = np.abs(np.gradient(np.gradient(mag_db))) + \
        np.abs(np.gradient(np.gradient(ph_detrend)))
    return float(f[np.argmax(score)])


def _fit_notch(data, outdir, out_name="testC_full_design_fit.json"):
    """Khalil et al. 2012 asymmetric hanger-notch fit of the complex S21(f):

        S21 = a e^{i alpha} e^{-2 pi i (f-f0) tau}
              * [ 1 - (Qr/|Qc|) e^{i phi} / (1 + 2i Qr (f/fr - 1)) ]

    Qr is the loaded/total Q; phi the Fano/asymmetry phase; |Qc| = Qr / (dip
    depth in the ideal case). Reports fr, Qr, Qc alongside PALACE's values.
    """
    from scipy.optimize import least_squares
    fitpath = os.path.join(outdir, out_name)
    f = np.asarray(data["freqs"], float)
    s21 = np.asarray(data["S21"])
    if s21.size < 8 or np.any(~np.isfinite(s21)):
        json.dump({"error": "no usable S21"}, open(fitpath, "w"))
        return
    f0 = f.mean()
    mag = np.abs(s21)
    mag_db = 20 * np.log10(mag)
    ph = np.unwrap(np.angle(s21))
    ph_ex = ph - np.polyval(np.polyfit(f, ph, 1), f)      # excess phase
    ph_swing_deg = np.degrees(ph_ex.max() - ph_ex.min())
    if (mag_db.max() - mag_db.min() < 0.1) and ph_swing_deg < 15:
        json.dump({"error": "no notch resolved",
                   "S21_db_span": float(mag_db.max() - mag_db.min()),
                   "S21_excess_phase_swing_deg": float(ph_swing_deg),
                   "S21_db_median": float(np.median(mag_db)),
                   "f_min_hz": float(f.min()), "f_max_hz": float(f.max()),
                   "n_points": int(f.size), "df_hz": float(np.median(np.diff(f))),
                   "note": ("S21 flat in both amplitude and phase over this "
                            "window/resolution; the hanger feature "
                            "(width ~fr/Qc) was not captured.")},
                  open(fitpath, "w"), indent=2)
        return
    # feature centre = sharpest point of |S21| dB OR excess phase (phase-only
    # feature has no amplitude dip)
    curv = np.abs(np.gradient(np.gradient(mag_db))) + \
        np.abs(np.gradient(np.gradient(ph_ex)))
    fr0 = float(f[np.argmax(curv)])
    # rough width from the excess-phase extent
    thr = 0.5 * (ph_ex.max() - ph_ex.min())
    inside = f[np.abs(ph_ex - np.median(ph_ex)) > thr]
    width = (inside.max() - inside.min()) if inside.size > 1 else (f[1] - f[0]) * 5
    Qr0 = float(np.clip(fr0 / max(width, f[1] - f[0]), 50, 5e5))
    a0 = float(np.median(mag[[0, 1, -2, -1]]))

    def model(p):
        fr, Qr, Qc, phi, a, alpha, tau = p
        env = a * np.exp(1j * alpha) * np.exp(-2j * np.pi * (f - f0) * tau)
        return env * (1 - (Qr / Qc) * np.exp(1j * phi) / (1 + 2j * Qr * (f / fr - 1)))

    def resid(p):
        d = model(p) - s21
        return np.concatenate([d.real, d.imag])

    best = None
    for Qc0 in (5e3, 1.4e4, 5e4):
        for phi0 in (0.0, 0.5, -0.5):
            p0 = [fr0, Qr0, Qc0, phi0, a0,
                  float(np.angle(s21[0])), 0.0]
            try:
                sol = least_squares(resid, p0, method="trf", max_nfev=40000,
                                    x_scale="jac")
                if best is None or sol.cost < best.cost:
                    best = sol
            except Exception:
                pass
    if best is None:
        json.dump({"error": "fit failed"}, open(fitpath, "w"))
        return
    fr, Qr, Qc, phi, a, alpha, tau = best.x
    Qr, Qc = abs(Qr), abs(Qc)
    Qi = 1.0 / (1.0 / Qr - np.cos(phi) / Qc) if (1 / Qr - np.cos(phi) / Qc) > 0 else None
    out = {
        "fr_hz": float(fr), "Qr_loaded": float(Qr), "Qc": float(Qc),
        "Qi_internal": (float(Qi) if Qi else None),
        "phi_rad": float(phi), "fit_cost": float(best.cost),
        "n_points": int(f.size),
        "palace_reference": {"fr_hz": 8.101e9, "Qc": 14078, "Qr": 533},
        "model": "S21 = a e^{i alpha} e^{-2pi i (f-f0) tau} "
                 "(1 - (Qr/|Qc|) e^{i phi}/(1 + 2i Qr (f/fr-1)))",
    }
    json.dump(out, open(fitpath, "w"), indent=2)
    # also dump the fitted curve for plotting
    np.savetxt(os.path.join(outdir, out_name.replace("_fit.json", "_fitcurve.csv")),
               np.column_stack([f, np.abs(model(best.x)),
                                20 * np.log10(np.abs(model(best.x)))]),
               delimiter=",", header="freq_hz,S21_mag_fit,S21_db_fit", comments="")


TESTS = {"A": run_A, "B": run_B, "C": run_C}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tests", nargs="+", help="A B C or 'all'")
    ap.add_argument("--graphical", action="store_true", default=True)
    ap.add_argument("--no-solve", dest="solve", action="store_false", default=True)
    ap.add_argument("--quick", action="store_true", help="tiny mesh/sweep shakeout")
    args = ap.parse_args()

    global QUICK
    QUICK = args.quick
    sel = ["A", "B", "C"] if "all" in args.tests else [t.upper() for t in args.tests]
    outdir = HERE
    os.makedirs(outdir, exist_ok=True)
    allres = {}
    for t in sel:
        print(f"\n{'='*60}\n  TEST {t}\n{'='*60}", flush=True)
        try:
            allres[t] = TESTS[t](outdir, args.solve, args.graphical)
        except Exception as e:
            allres[t] = {"error": repr(e), "traceback": traceback.format_exc()}
            print(allres[t]["traceback"], flush=True)
        with open(os.path.join(outdir, "_run_results.json"), "w") as f:
            json.dump(allres, f, indent=2)
        print(json.dumps(allres[t], indent=2, default=str), flush=True)


if __name__ == "__main__":
    main()
