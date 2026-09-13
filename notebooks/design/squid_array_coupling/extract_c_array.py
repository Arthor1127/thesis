"""
Derive the flux-coupling constant c for the SQUID-ARRAY design from the
epr_result.json files written by run_amr.py. Pure post-processing - no
qiskit_metal / SQDMetal / Palace - so it re-runs in a second without
re-solving.

------------------------------------------------------------------------
Reformulation of ../squid_resonator_epr/EPR_C_EXTRACTION.md for an array
------------------------------------------------------------------------
The reference design has ONE SQUID loop: one arm carries a lumped inductor
L_probe, the other is shorted, and Palace reports one participation p_probe
(fraction of the mode's inductive energy in that port). Then

    c_m = (1/2) sqrt( p_probe * f_m / (2 * E_probe/h) ),
    E_probe/h = phi0^2 / (L_probe h) = 163.4 GHz / L_probe[nH]

Here the array has N loops, each with BOTH junction gaps carrying L_J (no arm
shorted - SquidArray always carves both). For loop k the resonator's current
drives a circulating current I_k through the loop; both junctions carry I_k,
so the loop's total inductive energy is  L_J I_k^2 = (1/2)(2 L_J) I_k^2, i.e.
the loop behaves as the reference's single arm with

    L_eff = 2 L_J     (both junctions in series around the loop)
    p_loop_k = p_{2k} + p_{2k+1}          (the loop's two junction ports)
    E_eff/h  = phi0^2 / (2 L_J h) = 163.4 GHz / (2 L_J[nH])

Applying the reference formula loop-by-loop:

    c_k = (1/2) sqrt( |p_loop_k| * f_m / (2 * E_eff/h) ),   sign = sign(p_loop_k)

The per-loop mutual inductances M_k add in series (the resonator flux threads
every loop), and c_k proportional to M_k, so the array's total coupling is the
signed sum

    c_array = sum_k c_k        (sign(p_loop_k) tracks the resonator current
                                profile: loops on opposite sides of a current
                                node contribute with opposite sign)

Also reported: sum_k |c_k| (incoherent upper bound) and the dominant single
loop.

L_J convergence (analog of EPR_C_EXTRACTION.md Sec 6.1's two-point L_probe
fit): with L_loop the loop's own geometric self-inductance,

    c_meas(L_J) = c_inf * (2 L_J) / (2 L_J + L_loop)

Linearise 1/c = 1/c_inf + (L_loop/c_inf) * 1/(2 L_J). Two L_J values give
c_inf (the geometric answer) and L_loop. If they agree to mesh noise,
L_loop << 2 L_J and c_meas ~= c_inf directly.
"""
import argparse
import json
import math

PHI0_SQ_OVER_H_GHZ_NH = 163.4   # phi0^2 / h, GHz*nH  (EPR_C_EXTRACTION.md Sec 5)
C_SANITY_MAX = 0.1


def E_eff_over_h_ghz(L_J_nH):
    """Energy scale of the per-loop effective inductor 2*L_J."""
    return PHI0_SQ_OVER_H_GHZ_NH / (2.0 * L_J_nH)


def c_k_from_p(p_loop_k, f_ghz, L_J_nH):
    """Signed per-loop coupling. |p| under the sqrt, sign carried through."""
    if p_loop_k == 0.0:
        return 0.0
    mag = 0.5 * math.sqrt(abs(p_loop_k) * f_ghz / (2.0 * E_eff_over_h_ghz(L_J_nH)))
    return math.copysign(mag, p_loop_k)


def per_loop_participations(res_mode, n_squids):
    """Return the N per-loop participations p_loop_k, from p_per_loop if the
    run wrote it, else summing junction pairs from p_per_junction."""
    if "p_per_loop" in res_mode and len(res_mode["p_per_loop"]) == n_squids:
        return list(res_mode["p_per_loop"])
    pj = res_mode["p_per_junction"]
    return [pj[2 * k] + pj[2 * k + 1] for k in range(n_squids)]


# Only loops whose |participation| is a real fraction of the peak carry
# signal; the rest sit at the ~1e-9 numerical floor with random sign and,
# summed over ~20 loops, can swamp the true coherent c. Keep loops above
# this fraction of the peak per-loop |p|. 0.08 cleanly isolates the two
# real contributors (array near + far end); a looser 0.03 let a
# third, borderline loop in whose random sign flipped the L_J trend.
SIGNAL_FRACTION_OF_PEAK = 0.08


def analyse_run(result):
    """One epr_result.json -> the array coupling for its resonator mode."""
    rm = result["resonator_mode"]
    n = result["n_squids"]
    L_J = result["l_j_nh"]
    f = rm["f_ghz"]
    p_loop = per_loop_participations(rm, n)
    c_k = [c_k_from_p(p, f, L_J) for p in p_loop]

    peak = max(abs(p) for p in p_loop) or 1.0
    signal = [i for i in range(n) if abs(p_loop[i]) >= SIGNAL_FRACTION_OF_PEAK * peak]
    c_array = sum(c_k[i] for i in signal)                  # coherent, signal loops only
    c_array_all = sum(c_k)                                 # every loop (noise included)
    c_incoh = sum(abs(c_k[i]) for i in signal)             # incoherent upper bound
    kdom = max(range(n), key=lambda i: abs(c_k[i]))
    return {
        "n_squids": n,
        "l_j_nh": L_J,
        "f_ghz": f,
        "Q": rm.get("Q"),
        "p_array_signed": rm.get("p_array_signed"),
        "p_array_abs_sum": rm.get("p_array_abs_sum"),
        "c_array": c_array,
        "c_array_all_loops": c_array_all,
        "c_array_incoherent": c_incoh,
        "signal_loops": signal,
        "signal_loop_p": [p_loop[i] for i in signal],
        "dominant_loop": {"index": kdom, "c_k": c_k[kdom], "p_loop": p_loop[kdom]},
        "c_per_loop": c_k,
        "p_per_loop": p_loop,
        "amr_converged": result.get("amr", {}).get("frequency_converged_lt_0.5pct"),
        "final_mesh_elements": (result.get("final_mesh", {}).get("MeshElements")
                                or result.get("mesh", {}).get("elements")),
    }


def fit_c_inf_L_loop(L_J_list, c_list):
    """c_meas(L_J) = c_inf * (2 L_J)/(2 L_J + L_loop).
    Linearise: y = a + b x, y = 1/c, x = 1/(2 L_J), a = 1/c_inf, b = L_loop/c_inf.
    Closed form for 2 points, least squares for more."""
    xs = [1.0 / (2.0 * L) for L in L_J_list]
    ys = [1.0 / c for c in c_list]
    n = len(xs)
    if n < 2:
        return None
    if n == 2:
        b = (ys[1] - ys[0]) / (xs[1] - xs[0])
        a = ys[0] - b * xs[0]
        resid = None
    else:
        mx, my = sum(xs) / n, sum(ys) / n
        den = sum((x - mx) ** 2 for x in xs)
        if den == 0:
            return None
        b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den
        a = my - b * mx
        resid = math.sqrt(sum((a + b * x - y) ** 2 for x, y in zip(xs, ys)) / n)
    c_inf = 1.0 / a
    L_loop = b * c_inf
    return {"c_inf": c_inf, "L_loop_nH": L_loop, "residual": resid}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--result", action="append", required=True,
                    help="epr_result.json from run_amr.py. Repeat for every "
                         "(L_J, n_squids) run.")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    runs = []
    for p in args.result:
        d = json.load(open(p))
        if d.get("status") != "success":
            print(f"  SKIP {p}: status={d.get('status')}")
            continue
        a = analyse_run(d)
        a["path"] = p
        runs.append(a)

    checks = {}
    for r in runs:
        if abs(r["c_array"]) > C_SANITY_MAX:
            checks[f"{r['path']}:c_sane"] = {
                "pass": False,
                "detail": f"c_array={r['c_array']:.4f} > {C_SANITY_MAX}"}
        if r["amr_converged"] is False:
            checks[f"{r['path']}:amr_converged"] = {
                "pass": False, "detail": "resonator frequency not converged <0.5% - c is provisional"}

    # --- L_J convergence, per n_squids group ---
    by_n = {}
    for r in runs:
        by_n.setdefault(r["n_squids"], []).append(r)
    lj_convergence = []
    for n, grp in sorted(by_n.items()):
        grp = sorted(grp, key=lambda r: r["l_j_nh"])
        L_list = [r["l_j_nh"] for r in grp]
        c_list = [r["c_array"] for r in grp]
        entry = {"n_squids": n, "l_j_nh": L_list, "c_array": c_list}
        if len(set(L_list)) >= 2 and all(c != 0 for c in c_list):
            spread = (max(map(abs, c_list)) - min(map(abs, c_list))) / max(map(abs, c_list))
            entry["relative_spread"] = spread
            entry["fit"] = fit_c_inf_L_loop(L_list, [abs(c) for c in c_list])
            checks[f"n{n}:LJ_convergence"] = {
                "pass": spread <= 0.05,
                "detail": (f"c_array spread across L_J = {spread:.3f} "
                           f"({'L_loop negligible, c ~= c_inf' if spread <= 0.05 else 'use the c_inf/L_loop fit'})")}
        lj_convergence.append(entry)

    # --- coupling vs number of SQUIDs (at fixed L_J) ---
    by_lj = {}
    for r in runs:
        by_lj.setdefault(r["l_j_nh"], []).append(r)
    n_squids_trend = []
    for L_J, grp in sorted(by_lj.items()):
        grp = sorted(grp, key=lambda r: r["n_squids"])
        if len(grp) < 2:
            continue
        ns = [r["n_squids"] for r in grp]
        cs = [r["c_array"] for r in grp]
        monotonic = all(abs(cs[i]) <= abs(cs[i + 1]) + 1e-12 for i in range(len(cs) - 1))
        n_squids_trend.append({"l_j_nh": L_J, "n_squids": ns, "c_array": cs,
                               "increases_with_n": monotonic})
        checks[f"Lj{L_J}:coupling_grows_with_n"] = {
            "pass": monotonic,
            "detail": f"c_array vs n_squids={list(zip(ns, [round(c, 6) for c in cs]))}"}

    derived = {
        "method": "EPR_C_EXTRACTION.md reformulated for a SQUID array (see module docstring)",
        "runs": [{k: v for k, v in r.items() if k != "c_per_loop"} for r in runs],
        "c_per_loop_by_run": {r["path"]: r["c_per_loop"] for r in runs},
        "lj_convergence": lj_convergence,
        "n_squids_trend": n_squids_trend,
        "checks": checks,
        "all_passed": all(c["pass"] for c in checks.values()) if checks else None,
    }
    print(json.dumps({k: v for k, v in derived.items() if k != "c_per_loop_by_run"}, indent=2))
    if args.out:
        json.dump(derived, open(args.out, "w"), indent=2)
        print(f"\nWritten to {args.out}")


if __name__ == "__main__":
    main()
