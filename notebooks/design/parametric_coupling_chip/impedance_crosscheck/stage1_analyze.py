"""STAGE 1 analysis: read the Palace CPW benchmark port-S.csv and check how
close the lumped-port result sits to the documented Z0 = 56.02 ohm, at Order 2
vs Order 4.

Port 1 is R = 56.02 (matched to the line), so a perfect result is |S11| -> 0
and Zin = R*(1+S11)/(1-S11) -> 56.02. Any residual reflection is the
lumped-port-on-CPW discretisation error the Palace docs describe.
"""
import glob
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
R_PORT = 56.02  # examples/cpw/cpw_lumped_uniform.json


def _read_portS(path):
    # palace port-S.csv: "f (GHz)", "|S[1][1]| (dB)", "arg(S[1][1]) (deg.)", ...
    with open(path) as f:
        header = [h.strip() for h in f.readline().split(",")]
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    if data.ndim == 1:
        data = data[None, :]
    col = {h: i for i, h in enumerate(header)}
    f_ghz = data[:, col["f (GHz)"]]
    mag = 10 ** (data[:, col["|S[1][1]| (dB)"]] / 20.0)
    ph = np.deg2rad(data[:, col["arg(S[1][1]) (deg.)"]])
    s11 = mag * np.exp(1j * ph)
    return f_ghz, s11


def main():
    print(f"reference (Visser / Palace docs): Z0 = {R_PORT} ohm\n")
    hits = sorted(glob.glob(os.path.join(HERE, "stage1_out", "cpw_src", "postpro_*", "port-S.csv")))
    if not hits:
        raise SystemExit("no stage1_out/**/port-S.csv -- run stage1_run_palace_cpw.sh first")
    for path in hits:
        order = os.path.basename(os.path.dirname(path)).replace("postpro_", "")
        f_ghz, s11 = _read_portS(path)
        zin = R_PORT * (1 + s11) / (1 - s11)
        s11_db = 20 * np.log10(np.abs(s11))
        print(f"--- {order}  ({path}) ---")
        for fi, s, z in zip(f_ghz, s11_db, zin):
            print(f"  f={fi:6.2f} GHz   |S11|={s:7.2f} dB   Zin={z.real:6.2f}{z.imag:+.2f}j ohm"
                  f"   |Zin|={abs(z):6.2f}  ({100*(abs(z)/R_PORT-1):+.1f}%)")
        print(f"  mean |Zin| = {np.mean(np.abs(zin)):.2f} ohm "
              f"({100*(np.mean(np.abs(zin))/R_PORT-1):+.1f}% vs 56.02), "
              f"mean |S11| = {np.mean(s11_db):.1f} dB\n")


if __name__ == "__main__":
    main()
