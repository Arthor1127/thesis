"""Fit the JJ-array plasma-mode (phonon) dispersion and extract the effective
capacitances, given the known per-junction inductance L_J.

Lumped model: a chain of N junctions (N+1 islands). Each junction is L_J in
parallel with C_J,eff; each island has C_g to ground. Linearising the junctions,
the *bulk* dispersion (independent of how the chain terminates) is

    omega(k)^2 = (2/L_J) (1 - cos k) / [ C_g + 2 C_J,eff (1 - cos k) ]

The boundary conditions fix how k is quantised and the mode shape:

  * PERIODIC (Born-von Karman): travelling waves,  k_m = 2 pi m / N.
  * OPEN / free-free (our case - both chain ends float): standing waves. The
    free-end circuit equation reduces to the Neumann condition a_{-1} = a_0,
    so the modes are island-centred cosines  a_j ~ cos(k (j + 1/2))  with

        k_m = (m + delta) * pi / (N + 1),     m = 1 .. N

    An ideal uniform chain has delta = 0. The two end islands are not identical
    to the bulk (each touches only one junction and carries extra edge
    capacitance to ground), which shifts the standing-wave phase, so `delta`
    (a mode-index / boundary offset) is left as a global fit parameter instead
    of being fixed at 0.

We use this OPEN formula, not the periodic k_m = 2 pi m / N. The field maps
(mode 1 has a charge antinode at *both* ends) confirm free-free rather than a
grounded / quarter-wave termination, i.e. delta should come out small.

L_J is fixed (1.2 nH). Free parameters: C_g, C_J,eff, delta - shared across all
N. Box / package modes (sum_j |p_j| ~ 0) are dropped before fitting.
"""

import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
L_J = 1.2e-9              # H, fixed (eigenmodes.py --L-j-nH 1.2)
C_J_SET = 2.0e-15        # F, the lumped junction C put on each port
NS = [10, 23, 40, 100]
ARRAY_MODE_MIN_ABS_SUM = 0.5   # sum|p_j| ~ 1 -> array mode; ~ 0 -> box mode


def load_modes(n):
    """Return (array-mode freqs [Hz], dropped box-mode freqs [Hz])."""
    d = json.load(open(f"{HERE}/eig/n{n}/modes.json"))
    f = np.array([m["f_ghz"] for m in d["modes"]]) * 1e9
    ps = np.array([m.get("p_array_abs_sum", 0.0) for m in d["modes"]])
    keep = ps > ARRAY_MODE_MIN_ABS_SUM
    return np.sort(f[keep]), np.sort(f[~keep])


def omega_model(k, C_g, C_Jeff):
    """Bulk plasmon dispersion (rad/s)."""
    return np.sqrt((2.0 / L_J) * (1 - np.cos(k)) /
                   (C_g + 2.0 * C_Jeff * (1 - np.cos(k))))


def k_open(m, n, delta, denom_is_Np1=True):
    """Open (free-free) standing-wave quantisation with a fitted index offset."""
    denom = (n + 1) if denom_is_Np1 else n
    return (m + delta) * np.pi / denom


def _stack(data):
    Ns, Ms, Fs = [], [], []
    for n, fs in data.items():
        Ns.append(np.full(len(fs), n))
        Ms.append(np.arange(1, len(fs) + 1))
        Fs.append(fs)
    return map(np.concatenate, (Ns, Ms, Fs))


def global_fit(data, fit_delta=True, denom_is_Np1=True, x0=(1.0, 9.0, 0.0)):
    Ns, Ms, Fs = _stack(data)
    w = 2 * np.pi * Fs

    def resid(x):
        Cg, CJ = x[0] * 1e-15, x[1] * 1e-15
        delta = x[2] if fit_delta else 0.0
        return (omega_model(k_open(Ms, Ns, delta, denom_is_Np1), Cg, CJ) - w) / w

    sol = least_squares(resid, x0=list(x0),
                        bounds=([1e-3, 1e-3, -1.5], [1e4, 1e4, 1.5]))
    Cg, CJ, delta = sol.x[0] * 1e-15, sol.x[1] * 1e-15, (sol.x[2] if fit_delta else 0.0)
    return dict(C_g=Cg, C_Jeff=CJ, delta=delta,
                rms=float(np.sqrt(np.mean(resid(sol.x) ** 2))))


def per_n_fit(n, fs, delta):
    ms = np.arange(1, len(fs) + 1)
    w = 2 * np.pi * fs

    def resid(x):
        return (omega_model(k_open(ms, n, delta), x[0] * 1e-15, x[1] * 1e-15) - w) / w

    sol = least_squares(resid, x0=[1.0, 9.0], bounds=([1e-3, 1e-3], [1e4, 1e4]))
    return (sol.x[0] * 1e-15, sol.x[1] * 1e-15,
            float(np.sqrt(np.mean(resid(sol.x) ** 2))))


def main():
    data, box = {}, {}
    print(f"{'N':>4} {'array':>6} {'box':>4}   f_array (GHz)                         "
          f"  f_box (GHz)")
    for n in NS:
        fs, fb = load_modes(n)
        data[n], box[n] = fs, fb
        a_str = ", ".join(f"{x/1e9:.2f}" for x in fs[:5]) + \
                (f" ... {fs[-1]/1e9:.2f}" if len(fs) > 5 else "")
        b_str = ", ".join(f"{x/1e9:.1f}" for x in fb) if len(fb) else "-"
        print(f"{n:>4} {len(fs):>6} {len(fb):>4}   {a_str:<40}  {b_str}")
    print("  (box modes: sum_j|p_j| < 0.5 - no inductive energy in the array; "
          "dropped)")

    g_ideal = global_fit(data, fit_delta=False)
    g_delta = global_fit(data, fit_delta=True)
    g_dN = global_fit(data, fit_delta=True, denom_is_Np1=False)

    print("\n================  GLOBAL FITS  ================")
    print(f"{'variant':>28}  {'C_g (fF)':>9}  {'C_J,eff (fF)':>12}  "
          f"{'delta':>7}  {'rms':>9}")
    for name, g in [("open, k=m·pi/(N+1)  (delta=0)", g_ideal),
                    ("open, k=(m+delta)·pi/(N+1)", g_delta),
                    ("open, k=(m+delta)·pi/N", g_dN)]:
        print(f"{name:>28}  {g['C_g']*1e15:>9.3f}  {g['C_Jeff']*1e15:>12.3f}  "
              f"{g['delta']:>7.3f}  {g['rms']:>9.2e}")

    best = g_delta
    print(f"\n  per-N fits (delta fixed at global {best['delta']:.3f}):")
    print(f"  {'N':>4}  {'C_g (fF)':>9}  {'C_J,eff (fF)':>12}  "
          f"{'sqrt(L/Cg) (Ohm)':>16}  {'rms':>9}")
    for n in NS:
        Cg, CJ, rms = per_n_fit(n, data[n], best["delta"])
        print(f"  {n:>4}  {Cg*1e15:>9.3f}  {CJ*1e15:>12.3f}  "
              f"{np.sqrt(L_J/Cg):>16.1f}  {rms:>9.2e}")

    d = best["delta"]
    k100 = k_open(np.arange(1, 5), 100, d)
    slope = np.sum(k100 * 2 * np.pi * data[100][:4]) / np.sum(k100 ** 2)
    C_g_ac = 1.0 / (L_J * slope ** 2)
    k_top = k_open(len(data[10]), 10, d)
    omc = 1 - np.cos(k_top)
    C_Jeff_be = ((2.0 / L_J) * omc / (2 * np.pi * data[10][-1]) ** 2 - C_g_ac) / (2 * omc)
    print("\n================  ASYMPTOTIC CROSS-CHECKS  ================")
    print(f"  C_g   from N=100 acoustic slope (k->0)  = {C_g_ac*1e15:.3f} fF")
    print(f"  C_J,eff from N=10 band edge (k~pi)       = {C_Jeff_be*1e15:.3f} fF")

    C_g, C_Jeff = best["C_g"], best["C_Jeff"]
    f_p = 1 / (2 * np.pi * np.sqrt(L_J * C_Jeff))
    f_edge = omega_model(np.pi, C_g, C_Jeff) / (2 * np.pi)
    print("\n================  BEST GLOBAL FIT  (open BC)  ================")
    print(f"  L_J (fixed)              = {L_J*1e9:.3f} nH")
    print(f"  C_g   (island -> ground) = {C_g*1e15:.3f} fF")
    print(f"  C_J,eff (per junction)   = {C_Jeff*1e15:.3f} fF   "
          f"(port C_J was {C_J_SET*1e15:.1f} fF; rest is inter-island geometry)")
    print(f"  boundary offset delta    = {d:.3f}   "
          f"(k_m = (m+delta)·pi/(N+1); ideal free-free: delta = 0)")
    print(f"  rms relative residual    = {best['rms']:.2e}")
    print(f"  array wave impedance  sqrt(L_J/C_g)        = {np.sqrt(L_J/C_g):.1f} Ohm")
    print(f"  junction plasma freq  1/2pi sqrt(L_J C_Jeff) = {f_p/1e9:.2f} GHz")
    print(f"  band-edge freq (k=pi)                      = {f_edge/1e9:.2f} GHz")

    # ---- plot ----
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(19, 5))
    colors = {n: c for n, c in zip(NS, plt.cm.viridis(np.linspace(0, 0.9, len(NS))))}

    for n in NS:
        fs = data[n]
        ms = np.arange(1, len(fs) + 1)
        ax0.plot(ms, fs / 1e9, "o", ms=5, color=colors[n], label=f"N={n}")
        mm = np.linspace(1, len(fs), 200)
        ax0.plot(mm, omega_model(k_open(mm, n, d), C_g, C_Jeff) / (2 * np.pi) / 1e9,
                 "-", lw=1.2, color=colors[n])
        for fb in box[n]:
            ax0.plot(len(fs) + 0.5, fb / 1e9, "x", color=colors[n], ms=7)
    ax0.set_xlabel("mode index  m   (x = dropped box modes)")
    ax0.set_ylabel("mode frequency (GHz)")
    ax0.set_title("frequency vs mode index\n(lines: global open-BC dispersion fit)")
    ax0.legend(fontsize=8)
    ax0.grid(alpha=0.3)

    kk = np.linspace(1e-3, np.pi, 400)
    ax1.plot(kk / np.pi, omega_model(kk, C_g, C_Jeff) / (2 * np.pi) / 1e9,
             "k-", lw=1.5, label="dispersion fit")
    for n in NS:
        fs = data[n]
        ms = np.arange(1, len(fs) + 1)
        ax1.plot(k_open(ms, n, d) / np.pi, fs / 1e9, "o", ms=5,
                 color=colors[n], label=f"N={n}")
    ax1.set_xlabel(r"$k/\pi = (m+\delta)/(N+1)$")
    ax1.set_ylabel("mode frequency (GHz)")
    ax1.set_title("plasma-mode dispersion (open BC)\n"
                  fr"$C_g$={C_g*1e15:.2f} fF, $C_{{J,eff}}$={C_Jeff*1e15:.2f} fF, "
                  fr"$\delta$={d:.3f}")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)

    for n in NS:
        fs = data[n]
        ms = np.arange(1, len(fs) + 1)
        k = k_open(ms, n, d)
        model = omega_model(k, C_g, C_Jeff) / (2 * np.pi)
        ax2.plot(k / np.pi, 100 * (fs - model) / model, "o-", ms=4,
                 color=colors[n], label=f"N={n}")
    ax2.axhline(0, color="k", lw=0.8)
    ax2.set_xlabel(r"$k/\pi$")
    ax2.set_ylabel("(sim - fit) / fit  (%)")
    ax2.set_title("best-fit residuals")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    os.makedirs(f"{HERE}/figures", exist_ok=True)
    out = f"{HERE}/figures/dispersion_fit.png"
    fig.savefig(out, dpi=120, bbox_inches="tight")
    print(f"\nplot -> {out}")


if __name__ == "__main__":
    main()
