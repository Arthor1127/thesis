"""Assemble SUMMARY.md and comparison figures from the per-test result JSON /
CSV files written by run_ansys_tests.py. Safe to re-run after any subset of the
tests has completed."""
import glob
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

ANALYTIC = {
    "DesignHelper (conformal map)": 49.96,
    "SQDMetal CpwParams (grounded backplane)": 50.29,
    "PALACE capacitance-per-length": 47.55,
    "PALACE driven S-param |Zin|": 40.48,
}


def _load(name):
    p = os.path.join(HERE, f"{name}_result.json")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def _csv(name):
    p = os.path.join(HERE, f"{name}_sparams.csv")
    if not os.path.exists(p):
        return None
    return np.genfromtxt(p, delimiter=",", names=True)


def _zo(name):
    p = os.path.join(HERE, f"{name}_Zo.csv")
    if not os.path.exists(p):
        return None
    return np.genfromtxt(p, delimiter=",", names=True)


def fig_testA(a):
    z = _zo("testA_line_only")
    s = _csv("testA_line_only")
    if z is None:
        return
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    f = z["freq_hz"] / 1e9
    ax[0].plot(f, z["re_Zpi_port1"], "o-", label="Re Zpi (HFSS wave-port field)")
    if "re_Zin_fromS11_50" in z.dtype.names:
        ax[0].plot(f, z["re_Zin_fromS11_50"], "s-", label="Re Zin = 50(1+S11)/(1-S11)")
        ax[0].plot(f, z["im_Zin_fromS11_50"], "s--", alpha=.5, label="Im Zin")
    for v, lab in [(49.96, "analytic 49.96"), (50.29, "CpwParams 50.29"),
                   (40.48, "PALACE Zin 40.5")]:
        ax[0].axhline(v, ls=":", lw=1)
    ax[0].set_xlabel("freq (GHz)"); ax[0].set_ylabel("Z (ohm)")
    ax[0].set_title("Test A: bare CPW line impedance"); ax[0].legend(fontsize=8)
    ax[0].grid(alpha=.3)
    if s is not None:
        ax[1].plot(s["freq_hz"] / 1e9, s["S11_mag_db_50ohm"], "o-", label="S11 (50 Ohm renorm)")
        ax[1].plot(s["freq_hz"] / 1e9, s["S21_mag_db_50ohm"], "s-", label="S21 (50 Ohm renorm)")
        ax[1].set_xlabel("freq (GHz)"); ax[1].set_ylabel("dB")
        ax[1].set_title("Test A: S-parameters"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "testA_line_only.png"), dpi=130)
    plt.close(fig)


def fig_compare(a, b):
    fig, ax = plt.subplots(figsize=(8, 4))
    labels, vals, colors = [], [], []
    for k, v in ANALYTIC.items():
        labels.append(k); vals.append(v)
        colors.append("#888" if "PALACE driven" not in k else "#c0392b")
    if a:
        labels.append("HFSS Zpi (Test A)"); vals.append(a.get("Zpi_port1_mean_real") or 0)
        colors.append("#2471a3")
        labels.append("HFSS Zin/S11 (Test A)")
        vals.append(a.get("Zin_fromS11_50_abs_mean") or 0); colors.append("#5dade2")
    if b:
        labels.append("HFSS Zin/S11 (Test B, tapers)")
        vals.append(b.get("Zin_fromS11_50_abs_mean") or 0); colors.append("#a9cce3")
    c = _load("testC_full_design")
    if c:
        labels.append("HFSS Zin/S11 (Test C, full design)")
        vals.append(c.get("Zin_fromS11_50_abs_mean") or 0); colors.append("#d4e6f1")
    labels.append("PALACE Zin (Test B, tapers)"); vals.append(29.7); colors.append("#e6a0a0")
    labels.append("PALACE Zin (Test C, full design)"); vals.append(21.6); colors.append("#e6a0a0")
    y = np.arange(len(labels))
    ax.barh(y, vals, color=colors)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(50, color="k", ls="--", lw=1)
    for i, v in enumerate(vals):
        ax.text(v + 0.3, i, f"{v:.1f}", va="center", fontsize=8)
    ax.set_xlabel("impedance (ohm)")
    ax.set_title("CPW characteristic impedance / Zin -- methods compared")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "comparison.png"), dpi=130)
    plt.close(fig)


def fig_testC(c):
    s = _csv("testC_full_design")
    fitp = os.path.join(HERE, "testC_full_design_fit.json")
    if s is None:
        return
    fig, ax = plt.subplots(figsize=(8, 4))
    y = s["S21_mag_db_50ohm"]
    ax.plot(s["freq_hz"] / 1e9, y, "-", lw=1, label="driven |S21|")
    ax.set_xlabel("freq (GHz)"); ax.set_ylabel("|S21| (dB)")
    ax.set_title("Test C: driven |S21| (flat) + eigenmode f_r"); ax.grid(alpha=.3)
    span = y.max() - y.min()
    ax.set_ylim(y.min() - 3, y.max() + 3)
    eig = os.path.join(HERE, "testC_eigenmodes.json")
    if os.path.exists(eig):
        with open(eig) as f:
            ej = json.load(f)
        ms = ej.get("modes", [])
        if ms and ms[0].get("f_hz"):
            fr = ms[0]["f_hz"] / 1e9
            ax.axvline(fr, color="r", ls="--",
                       label=f"eigenmode f_r = {fr:.3f} GHz (Q_i≈{ms[0].get('Q',0):.0f})")
    ax.set_ylim(-2, 1)
    ax.text(0.5, 0.15, "driven S21 essentially flat (raw |S11| = 0.054 +/- 0.002,\n"
            "no resonance) — coupling too weak for the driven mesh to resolve;\n"
            "f_r is from the eigenmode solve. The 0.7 dB blip at 8.43 GHz is a\n"
            "port-Zo renormalisation artefact near the mode.",
            transform=ax.transAxes, ha="center", va="center", fontsize=8, color="#666")
    ax.legend(fontsize=8, loc="upper right")
    fcv = os.path.join(HERE, "testC_full_design_fitcurve.csv")
    if os.path.exists(fcv):
        fc = np.genfromtxt(fcv, delimiter=",", names=True)
        ax.plot(fc["freq_hz"] / 1e9, fc["S21_db_fit"], "r--", lw=1, label="Khalil fit")
    if os.path.exists(fitp):
        with open(fitp) as f:
            fj = json.load(f)
        if "fr_hz" in fj:
            qr = fj.get("Qr_loaded", fj.get("Qr"))
            ax.axvline(fj["fr_hz"] / 1e9, color="r", ls=":",
                       label=f"fr={fj['fr_hz']/1e9:.4f} GHz Qc={fj['Qc']:.0f} Qr={qr:.0f}")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "testC_full_design.png"), dpi=130)
    plt.close(fig)


def main():
    a = _load("testA_line_only")
    b = _load("testB_launchpads")
    c = _load("testC_full_design")
    fig_testA(a)
    fig_compare(a, b)
    fig_testC(c)

    L = []
    w = L.append
    w("# ANSYS HFSS cross-check of the CPW impedance discrepancy -- results\n")
    w("See `../HANDOFF.md` for the question. All runs: ANSYS Electronics Desktop "
      "2024 R2, HFSS Driven-Modal, **wave ports** with an integration line "
      "across the CPW gap (characteristic impedance `Zo`/`Zpi` taken directly "
      "from the 2D port field solution), zero-thickness Perfect-E metal, "
      "silicon `er=11.45` lossless, grounded Perfect-E backplane 500 um below "
      "the trace. S-parameters exported non-renormalized then renormalized to a "
      "flat 50 Ohm in post so `Zin = 50 (1+S11)/(1-S11)` matches PALACE's "
      "convention.\n")

    def block(tag, r, extra=""):
        if not r:
            w(f"## {tag}\n\n_(not yet run)_\n")
            return
        w(f"## {tag}\n")
        zpi = r.get("Zpi_port1_mean_real")
        zpi_i = r.get("Zpi_port1_mean_imag") or 0.0
        zpi_s = r.get("Zpi_port1_std_real") or 0.0
        if zpi is not None:
            w(f"- **HFSS field Zpi (port 1): {zpi:.2f} {zpi_i:+.2f}j Ohm** "
              f"(std over band {zpi_s:.3f})")
        else:
            w("- HFSS field Zpi: unavailable")
        if r.get("eps_eff_port1_mean"):
            w(f"- port-mode eps_eff: {r['eps_eff_port1_mean']:.3f} "
              f"(analytic ~5.49)")
        w(f"- Zin from S11 (renorm 50 Ohm), band mean: "
          f"{r['Zin_fromS11_50_mean_real']:.2f} {r['Zin_fromS11_50_mean_imag']:+.2f}j Ohm "
          f"(|Zin| {r['Zin_fromS11_50_abs_mean']:.2f})")
        w(f"- mean |S11| (50 Ohm): {r['S11_50_mean_db']:.1f} dB;  "
          f"non-renorm |S11| (ref = Zo): {r.get('S11_raw_mean_abs', float('nan')):.3f}")
        cv = r.get("convergence", {}).get("raw")
        if cv:
            w(f"- convergence tail: `{cv[-1]}`")
        if r.get("warnings"):
            w(f"- warnings: {r['warnings']}")
        w(extra)
        w("")

    block("Test A -- bare CPW line", a)
    block("Test B -- CPW + launchpad tapers", b)
    block("Test C -- full hanger-resonator design", c)

    fitp = os.path.join(HERE, "testC_full_design_fit.json")
    if os.path.exists(fitp):
        with open(fitp) as f:
            fj = json.load(f)
        w("### Test C notch fit (Khalil et al. 2012 asymmetric hanger)\n")
        w(f"```\n{json.dumps(fj, indent=2)}\n```\n")

    w("## Deviations from the handoff\n")
    w("- **Compact domain box** (a few hundred um of silicon/vacuum/ground "
      "around the structure) rather than the full 4.8 x 2 mm chip. Justified "
      "by the handoff's own finding that lateral ground width (0.1x-3x) and "
      "vacuum height (1x-12x) have zero effect; a compact uniform cross-section "
      "is the standard choice for wave-port field accuracy.\n")
    w("- **No wirebonds / airbridges** modeled (planar only). The CPW mode is "
      "excited symmetrically by the wave port so the slotline mode the "
      "wirebonds would suppress is not excited.\n")
    w("- **Zero-thickness Perfect-E metal** (as requested, matching SQDMetal's "
      "quasi-3D treatment).\n")
    w("- Wave ports sit on the domain end faces; the modeled ground therefore "
      "ends at the port plane instead of extending past it.\n")

    with open(os.path.join(HERE, "_autosummary.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L))
    print("wrote figures + _autosummary.md (SUMMARY.md is the authoritative one)")


if __name__ == "__main__":
    main()
