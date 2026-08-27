"""Compare isolated-wire vs CPW center-conductor inductance as a function
of trace width, to check whether widening the SQUID arm reduces its
self-inductance. Not wired into any design - standalone check."""
import numpy as np
import pint
from scipy.special import ellipk

ureg = pint.UnitRegistry()
mu_0 = ureg.Quantity(1, ureg.mu_0).to('H/m')


def line_inductance_isolated(l, w, t):
    l = l.to("um"); w = w.to("um"); t = t.to("um")
    arg_1 = (2.0 * l / (w + t)).magnitude
    arg_2 = 0.2235 * ((w + t) / l).magnitude
    L = mu_0 * l * (np.log(arg_1) + 0.5 + arg_2) / (2.0 * np.pi)
    return L.to("henry")


def line_inductance_cpw(l, w, g):
    w = w.to("um").magnitude
    g = g.to("um").magnitude
    k = w / (w + 2.0 * g)
    kp = np.sqrt(1.0 - k**2)
    Kk = ellipk(k**2)
    Kkp = ellipk(kp**2)
    L_prime = mu_0 * 0.25 * (Kkp / Kk)
    return (L_prime * l.to("um")).to("henry")


if __name__ == "__main__":
    l = ureg.Quantity(1200, "um")   # arm length, l1+l2
    t = ureg.Quantity(100, "nm")    # film thickness
    g = ureg.Quantity(4, "um")      # gap to ground (matches squid_options.g1)

    widths_um = np.linspace(2, 30, 50)
    L_isolated = [line_inductance_isolated(l, ureg.Quantity(w, "um"), t).to("pH").magnitude
                  for w in widths_um]
    L_cpw = [line_inductance_cpw(l, ureg.Quantity(w, "um"), g).to("pH").magnitude
             for w in widths_um]

    print(f"{'w (um)':>8} {'L_isolated (pH)':>16} {'L_cpw (pH)':>12}")
    for w, li, lc in zip(widths_um[::5], L_isolated[::5], L_cpw[::5]):
        print(f"{w:>8.1f} {li:>16.2f} {lc:>12.2f}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(widths_um, L_isolated, label="isolated wire")
        ax.plot(widths_um, L_cpw, label="CPW (g=4um)")
        ax.set_xlabel("arm width w (um)")
        ax.set_ylabel("inductance (pH)")
        ax.set_title(f"Arm self-inductance vs width (l={l})")
        ax.legend()
        fig.savefig("width_vs_inductance.png", dpi=150)
        print("Saved width_vs_inductance.png")
    except ImportError:
        pass
