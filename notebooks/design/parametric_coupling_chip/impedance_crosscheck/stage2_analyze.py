"""STAGE 2 analysis: exact per-frequency recovery of Z0 and gamma for the Si
CPW through line from its (S11, S21), compared to the analytic Visser Z0.

Uniform-line (Nicolson-Ross-Weir style) inversion, reference impedance R:

    K   = (S11^2 - S21^2 + 1) / (2 S11)
    Gam = K +/- sqrt(K^2 - 1)              # pick |Gam| <= 1
    Z0  = R (1 + Gam) / (1 - Gam)
    T   = (S11 + S21 - Gam) / (1 - (S11 + S21) Gam)   # = exp(-gamma * l)
    gamma = -ln(T) / l  ->  eps_eff = (c * Im(gamma) / omega)^2

Exact for a matched-reference uniform line -- no fitting. When Palace gets the
CPW right, S11 -> 0, the line is matched, and Z0 -> R trivially; a real
mismatch (the ~40 ohm SQDMetal saw => |S11| ~ -19 dB) is recovered directly.
"""
import glob
import os
import numpy as np

from analytic_cpw import cpw_impedance

HERE = os.path.dirname(os.path.abspath(__file__))
R_PORT = 50.0
TRACE_UM = float(os.environ.get("TRACE_UM", "20.0"))
GAP_UM = float(os.environ.get("GAP_UM", "11.5"))
C0 = 299792458.0
LINE_LEN_M = 4000e-6      # keep in sync with stage2_run.sh LINE_LEN_UM


def _read_portS(path):
    with open(path) as f:
        header = [h.strip() for h in f.readline().split(",")]
    d = np.genfromtxt(path, delimiter=",", skip_header=1)
    if d.ndim == 1:
        d = d[None, :]
    col = {h: i for i, h in enumerate(header)}

    def cplx(which):
        mag = 10 ** (d[:, col[f"|S{which}| (dB)"]] / 20.0)
        ph = np.deg2rad(d[:, col[f"arg(S{which}) (deg.)"]])
        return mag * np.exp(1j * ph)

    return d[:, col["f (GHz)"]] * 1e9, cplx("[1][1]"), cplx("[2][1]")


def invert_line(f_hz, s11, s21, length_m):
    K = (s11**2 - s21**2 + 1.0) / (2.0 * s11)
    root = np.sqrt(K**2 - 1.0)
    g1, g2 = K + root, K - root
    gam = np.where(np.abs(g1) <= 1.0, g1, g2)
    z0 = R_PORT * (1.0 + gam) / (1.0 - gam)

    t = (s11 + s21 - gam) / (1.0 - (s11 + s21) * gam)
    # gamma*l = -ln(T); Im must be unwrapped in frequency for eps_eff
    gl = -np.log(t)
    beta_l = np.unwrap(gl.imag)
    beta = beta_l / length_m
    omega = 2 * np.pi * f_hz
    eps_eff = (C0 * beta / omega) ** 2
    alpha_np_m = gl.real / length_m
    return z0, eps_eff, alpha_np_m, gam


def main():
    Z0_an, eps_an = cpw_impedance(TRACE_UM, GAP_UM, 525.0, 11.45)
    print(f"analytic (Visser)  Si w={TRACE_UM:g} s={GAP_UM:g} h=525 epsr=11.45 -> "
          f"Z0 = {Z0_an:.2f} ohm   eps_eff = {eps_an:.3f}\n")

    outdir = os.environ.get("OUTDIR", "stage2_out")
    hits = sorted(glob.glob(os.path.join(HERE, outdir, "postpro_o*", "port-S.csv")))
    if not hits:
        raise SystemExit(f"no {outdir}/postpro_o*/port-S.csv -- run stage2_run.sh first")

    for path in hits:
        order = os.path.basename(os.path.dirname(path)).replace("postpro_o", "")
        f_hz, s11, s21 = _read_portS(path)
        z0, eps_eff, alpha, gam = invert_line(f_hz, s11, s21, LINE_LEN_M)

        s11_db = 20 * np.log10(np.abs(s11))
        # drop points where |S11| is at the numerical floor (inversion ill-posed there)
        ok = np.abs(s11) > 3e-3
        z0_ok = z0[ok].real
        print(f"--- Order {order}  ({ok.sum()}/{len(f_hz)} usable freqs) ---")
        for fi, sd, z, e in zip(f_hz, s11_db, z0, eps_eff):
            flag = "" if abs(z) > 0 else ""
            print(f"  f={fi/1e9:5.1f} GHz  |S11|={sd:6.1f} dB  Z0={z.real:6.2f}{z.imag:+5.2f}j  eps_eff={e:6.3f}")
        if ok.sum():
            print(f"  --> Z0  = {np.mean(z0_ok):6.2f} +- {np.std(z0_ok):.2f} ohm   "
                  f"({100*(np.mean(z0_ok)/Z0_an-1):+.1f}% vs analytic {Z0_an:.1f})")
            print(f"      eps_eff = {np.mean(eps_eff[ok]):.3f} +- {np.std(eps_eff[ok]):.3f}   "
                  f"({100*(np.mean(eps_eff[ok])/eps_an-1):+.1f}% vs analytic {eps_an:.3f})")
        else:
            print("  all |S11| below 3e-3 -> line is well matched; Z0 ~ R = 50 ohm to within the floor")
        print(f"      |S11| envelope: {s11_db.min():.1f} .. {s11_db.max():.1f} dB\n")


if __name__ == "__main__":
    main()
