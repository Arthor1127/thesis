"""Analytic characteristic impedance of a coplanar waveguide.

Direct port of `cpw_impedance(; w, s, h, epsr)` from the Palace source example
`examples/cpw/cpw.jl` (H. J. Visser, *Antenna Theory and Applications*, Wiley
2012, p. 259) -- verified to match the Julia to 3 d.p.

Palace's CPW benchmark quotes Z0 = 56.02 Ohm; that value comes from the
IN-PLANE sapphire permittivity epsr = 9.3 (perpendicular to the C-axis), not
the 11.5 surface-normal value -- the CPW fields are mostly in-plane. With
epsr = 9.3 this function returns 56.02; the ports in cpw_lumped_uniform.json
are set to R = 56.02 to match.

  w   trace width           [um]
  s   gap (trace<->ground)  [um]
  h   substrate thickness   [um]
  epsr relative permittivity, surface-normal direction
"""
import numpy as np


def _s(x):
    r = np.sqrt(1.0 + x)
    q = (4.0 * x) ** 0.25
    return np.log(2.0 * (r + q) / (r - q))


def _K_over_Kp(k, kp):
    if k >= 1.0 / np.sqrt(2.0):
        return _s(k) / (2.0 * np.pi)
    return 2.0 * np.pi / _s(kp)


def cpw_impedance(w=30.0, s=18.0, h=500.0, epsr=11.5):
    k = w / (w + 2.0 * s)
    k1 = np.sinh(np.pi * w / (4.0 * h)) / np.sinh(np.pi * (w + 2.0 * s) / (4.0 * h))
    kp = np.sqrt(1.0 - k**2)
    k1p = np.sqrt(1.0 - k1**2)
    KoKp = _K_over_Kp(k, kp)
    K1oK1p = _K_over_Kp(k1, k1p)
    eps_eff = 1.0 + ((epsr - 1.0) / 2.0) * K1oK1p / KoKp
    Z0 = 30.0 * np.pi / (KoKp * np.sqrt(eps_eff))
    return Z0, eps_eff


if __name__ == "__main__":
    # Palace benchmark cross-check -- epsr=9.3 must reproduce their 56.02:
    z, e = cpw_impedance(30.0, 18.0, 500.0, 9.3)
    print(f"Palace cpw.jl benchmark  w=30 s=18 h=500 epsr=9.3  -> Z0={z:.2f} ohm  eps_eff={e:.3f}  (docs: 56.02)")
    # This thesis' readout CPW cross-section (Si, isotropic):
    z, e = cpw_impedance(20.0, 11.5, 525.0, 11.45)
    print(f"thesis readout CPW       w=20 s=11.5 h=525 epsr=11.45 -> Z0={z:.2f} ohm  eps_eff={e:.3f}")
