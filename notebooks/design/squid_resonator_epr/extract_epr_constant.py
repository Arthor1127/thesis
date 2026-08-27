"""
Derive the flux-coupling constant c from one or more epr_result.json files
written by build_and_eigenmode.py. Pure post-processing - no qiskit_metal,
SQDMetal, or Palace dependency - so this can be re-run cheaply without
re-solving. See EPR_C_EXTRACTION.md for the full physics derivation.
"""
import argparse
import json
import math

PHI0_OVER_H_GHZ_NH = 163.4  # phi0^2 / h, in GHz*nH (EPR_C_EXTRACTION.md Sec 5)
C_SANITY_MAX = 0.1          # doc Sec 1/7: c > 0.1 means something is wrong


def E_probe_over_h_ghz(L_probe_nH):
    return PHI0_OVER_H_GHZ_NH / L_probe_nH


def c_from_participation(p_probe, f_ghz, L_probe_nH):
    """EPR_C_EXTRACTION.md Sec 5: c_m = 1/2 * sqrt(p_probe * f_m / (2 * E_probe/h))."""
    if p_probe <= 0:
        return 0.0
    return 0.5 * math.sqrt(p_probe * f_ghz / (2 * E_probe_over_h_ghz(L_probe_nH)))


def load_result(path):
    with open(path) as f:
        data = json.load(f)
    if data.get("status") != "success":
        raise ValueError(f"{path}: not a successful run (status={data.get('status')})")
    return data


def resonator_modes(result):
    return [m for m in result["modes"] if m["assignment"] != "transmon"]


def match_modes_across_runs(runs):
    """Match resonator modes across runs by nearest frequency (L_probe
    perturbs eigenfrequencies slightly, so ordinal position isn't safe -
    EPR_C_EXTRACTION.md Sec 6.1). Returns a list of per-mode-group records,
    one per mode in the first run, each holding the matched mode from every
    run (or None if no reasonable match was found)."""
    base_modes = resonator_modes(runs[0])
    groups = []
    for base_mode in base_modes:
        group = {0: base_mode}
        for run_idx in range(1, len(runs)):
            candidates = resonator_modes(runs[run_idx])
            if not candidates:
                group[run_idx] = None
                continue
            closest = min(candidates, key=lambda m: abs(m["f_ghz"] - base_mode["f_ghz"]))
            group[run_idx] = closest
        groups.append(group)
    return groups


def fit_c_inf_L_loop(L_probe_list, c_meas_list):
    """Solve c_meas(L_probe) = c_inf * L_probe / (L_probe + L_loop).
    Closed-form for exactly 2 points; least-squares (via linearization) for
    more, per EPR_C_EXTRACTION.md Sec 6.1.

    Linearize: 1/c_meas = (1/c_inf) + (L_loop/c_inf) * (1/L_probe)
    i.e. y = a + b*x, with y=1/c_meas, x=1/L_probe, a=1/c_inf, b=L_loop/c_inf.
    """
    xs = [1.0 / L for L in L_probe_list]
    ys = [1.0 / c for c in c_meas_list]
    n = len(xs)
    if n < 2:
        return None
    if n == 2:
        x1, x2 = xs
        y1, y2 = ys
        b = (y2 - y1) / (x2 - x1)
        a = y1 - b * x1
    else:
        mean_x = sum(xs) / n
        mean_y = sum(ys) / n
        num = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
        den = sum((x - mean_x) ** 2 for x in xs)
        b = num / den
        a = mean_y - b * mean_x
    c_inf = 1.0 / a
    L_loop = b * c_inf
    residual = None
    if n > 2:
        residual = math.sqrt(sum((a + b * x - y) ** 2 for x, y in zip(xs, ys)) / n)
    return {"c_inf": c_inf, "L_loop_nH": L_loop, "residual": residual}


def mesh_convergence_check(mesh_sweep_runs, rel_tol=0.05):
    """EPR_C_EXTRACTION.md Sec 6.2: report p_probe vs mesh element count,
    flag if not plateaued between the two finest (largest element count)
    points."""
    points = []
    for run in mesh_sweep_runs:
        modes = resonator_modes(run)
        if not modes:
            continue
        # use the mode with the largest p_probe as the one of interest
        mode = max(modes, key=lambda m: m["p"]["probe"])
        points.append({
            "elements": run["mesh"]["elements"],
            "junction_min_size_um": run["mesh"]["junction_min_size_um"],
            "p_probe": mode["p"]["probe"],
            "f_ghz": mode["f_ghz"],
        })
    points.sort(key=lambda pt: pt["elements"])
    plateaued = None
    rel_change = None
    if len(points) >= 2:
        p_coarse, p_fine = points[-2]["p_probe"], points[-1]["p_probe"]
        rel_change = abs(p_fine - p_coarse) / p_fine if p_fine else float("inf")
        plateaued = rel_change <= rel_tol
    return {"points": points, "rel_change_finest_two": rel_change, "plateaued": plateaued}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--result", action="append", required=True,
                     help="Path to an epr_result.json. Repeat for the mandatory "
                          "Sec 6.1 two-point L_probe fit (>=2 required for a fit).")
    ap.add_argument("--mesh-sweep-result", action="append", default=[],
                     help="Path to an epr_result.json from a junction-mesh-density "
                          "sweep at fixed L_probe, for the Sec 6.2 convergence check "
                          "(separate from --result, which sweeps L_probe).")
    ap.add_argument("--out", default=None, help="Write the derived JSON here too.")
    args = ap.parse_args()

    runs = [load_result(p) for p in args.result]
    L_probe_list = [r["l_probe_nh"] for r in runs]

    checks = {}
    per_run_modes = []
    for run, path in zip(runs, args.result):
        modes_out = []
        for mode in resonator_modes(run):
            c_m = c_from_participation(mode["p"]["probe"], mode["f_ghz"], run["l_probe_nh"])
            modes_out.append({**mode, "c_m": c_m})
            if mode["p"]["probe"] <= 0:
                checks[f"{path}:m{mode['m']}:p_probe_nonzero"] = {
                    "pass": False, "detail": "p_probe == 0, no circulating current path"}
            if c_m > C_SANITY_MAX:
                checks[f"{path}:m{mode['m']}:c_sane"] = {
                    "pass": False, "detail": f"c={c_m:.4f} > {C_SANITY_MAX} - implausible, "
                                              "check bridging/port setup"}
        per_run_modes.append({"path": path, "l_probe_nh": run["l_probe_nh"], "modes": modes_out})

    fit_result = None
    if len(runs) >= 2:
        groups = match_modes_across_runs(runs)
        fit_result = []
        for group in groups:
            base = group[0]
            L_list, c_list = [], []
            for run_idx, run in enumerate(runs):
                matched = group[run_idx]
                if matched is None:
                    continue
                c_m = c_from_participation(matched["p"]["probe"], matched["f_ghz"],
                                            run["l_probe_nh"])
                L_list.append(run["l_probe_nh"])
                c_list.append(c_m)
            entry = {"matched_base_f_ghz": base["f_ghz"], "L_probe_nh": L_list, "c_meas": c_list}
            if len(L_list) >= 2 and all(c > 0 for c in c_list):
                fit = fit_c_inf_L_loop(L_list, c_list)
                entry.update(fit)
                spread = (max(c_list) - min(c_list)) / max(c_list)
                entry["relative_spread"] = spread
                key = f"m~{base['f_ghz']:.3f}GHz:loop_negligible"
                checks[key] = {
                    "pass": spread <= 0.05,
                    "detail": (f"c_meas relative spread across L_probe points = {spread:.3f} "
                               f"({'loop inductance negligible' if spread <= 0.05 else 'loop is screening, use c_inf/L_loop fit'})")
                }
                if fit["L_loop_nH"] is not None and fit["L_loop_nH"] > 0.1:
                    checks[key.replace("loop_negligible", "L_loop_magnitude")] = {
                        "pass": False,
                        "detail": (f"L_loop={fit['L_loop_nH']*1000:.1f}pH is in the hundreds-of-pH "
                                   "range - Sec 3/8 escalation trigger (missing kinetic inductance "
                                   "may matter here)")
                    }
            fit_result.append(entry)

    mesh_check = None
    if args.mesh_sweep_result:
        mesh_runs = [load_result(p) for p in args.mesh_sweep_result]
        mesh_check = mesh_convergence_check(mesh_runs)
        if mesh_check["plateaued"] is not None:
            checks["mesh_convergence"] = {
                "pass": mesh_check["plateaued"],
                "detail": (f"p_probe relative change between finest two mesh points = "
                           f"{mesh_check['rel_change_finest_two']:.3f}")
            }

    derived = {
        "runs": per_run_modes,
        "two_point_fit": fit_result,
        "mesh_convergence": mesh_check,
        "checks": checks,
        "all_passed": all(c["pass"] for c in checks.values()) if checks else None,
    }

    print(json.dumps(derived, indent=2))
    if args.out:
        with open(args.out, "w") as f:
            json.dump(derived, f, indent=2)
        print(f"\nWritten to {args.out}")


if __name__ == "__main__":
    main()
