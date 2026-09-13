"""Recompute the 50-ohm-referenced S-parameters and Zin for the Test A runs
from the raw (Zpi-referenced) Touchstone + the Zo CSV, without re-solving.

HFSS exports non-renormalized S-parameters (reference = each port's Zpi(f)),
but the Touchstone header says 50 ohm, so scikit-rf's renormalize() is a no-op
unless we first set the true reference impedance. This rebuilds:
    <name>_50ohm.s2p        S renormalized to a flat 50 ohm
    <name>_Zo.csv           freq, Zpi(port1), Zin=50*(1+S11_50)/(1-S11_50)
and prints the band-mean numbers.
"""
import glob
import json
import os

import numpy as np
import skrf

HERE = os.path.dirname(os.path.abspath(__file__))


def reprocess(tag):
    s2p = os.path.join(HERE, f"{tag}.s2p")
    zocsv = os.path.join(HERE, f"{tag}_Zo.csv")
    if not (os.path.exists(s2p) and os.path.exists(zocsv)):
        return None
    z = np.genfromtxt(zocsv, delimiter=",", names=True)
    f = z["freq_hz"]
    zpi = z["re_Zpi_port1"] + 1j * z["im_Zpi_port1"]

    ntwk = skrf.Network(s2p)
    assert np.allclose(ntwk.f, f, rtol=1e-6), (tag, "freq mismatch")
    S11_raw = ntwk.s[:, 0, 0]
    # true reference = Zpi(f) on both ports (symmetric uniform line: port2 ~ port1)
    ntwk.z0 = np.column_stack([zpi, zpi])
    ntwk.renormalize(50.0)
    S11 = ntwk.s[:, 0, 0]
    S21 = ntwk.s[:, 1, 0]
    Zin = 50.0 * (1 + S11) / (1 - S11)

    ntwk.write_touchstone(os.path.join(HERE, f"{tag}_50ohm.s2p"), write_z0=True)
    np.savetxt(zocsv,
               np.column_stack([f, zpi.real, zpi.imag, Zin.real, Zin.imag,
                                np.abs(S11_raw), np.abs(S11)]),
               delimiter=",",
               header=("freq_hz,re_Zpi_port1,im_Zpi_port1,"
                       "re_Zin_fromS11_50,im_Zin_fromS11_50,"
                       "abs_S11_raw_Zpiref,abs_S11_50ohm"), comments="")

    m = (f >= 6e9) & (f <= 10e9)
    out = {
        "tag": tag,
        "Zpi_port1_mean": float(np.mean(zpi.real[m])),
        "Zin_fromS11_50_mean_real": float(np.mean(Zin.real[m])),
        "Zin_fromS11_50_mean_imag": float(np.mean(Zin.imag[m])),
        "Zin_fromS11_50_abs_mean": float(np.mean(np.abs(Zin[m]))),
        "abs_S11_50ohm_mean": float(np.mean(np.abs(S11[m]))),
        "abs_S11_raw_mean": float(np.mean(np.abs(S11_raw[m]))),
    }
    rj = os.path.join(HERE, f"{tag}_result.json")
    if os.path.exists(rj):
        d = json.load(open(rj))
        d.update({k: out[k] for k in out if k != "tag"})
        json.dump(d, open(rj, "w"), indent=2)
    return out


if __name__ == "__main__":
    tags = [os.path.basename(p)[:-4] for p in sorted(glob.glob(os.path.join(HERE, "testA*.s2p")))
            if not p.endswith("_50ohm.s2p")]
    for t in tags:
        r = reprocess(t)
        if r:
            print(f"{t:24s} Zpi={r['Zpi_port1_mean']:.2f}  "
                  f"Zin(50)={r['Zin_fromS11_50_mean_real']:.2f}{r['Zin_fromS11_50_mean_imag']:+.2f}j  "
                  f"|Zin(50)|={r['Zin_fromS11_50_abs_mean']:.2f}  "
                  f"|S11@50|={r['abs_S11_50ohm_mean']:.4f}  |S11@Zpi|={r['abs_S11_raw_mean']:.4f}")
