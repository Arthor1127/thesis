"""Write the sweep-grid JSONs.

Each grid fixes two of the four knobs  (omega_C, chi_C, omega_r, gamma)  as the
colormap axes (y vertical, x horizontal) and sweeps the other two discretely,
one sweep per (disc1, disc2) pair -> 8 x 8 = 64 sweeps.

  grid_chiC.json      axes (omega_C vert, chi_C   horiz)   disc: omega_r x gamma
  grid_omegar.json    axes (omega_C vert, omega_r horiz)   disc: chi_C   x gamma
  grid_gammaomegar.json  axes (gamma vert, omega_r horiz)  disc: omega_C x chi_C
"""
import json
import numpy as np

KNOBS = ("omega_C", "chi_C", "omega_r", "gamma")

BASE = dict(
    Omega_J=1.0, Chi_J=0.0,
    d_ph=4, Nk=96, s_A=-0.25, s_B=0.25,
    gamma_values=list(np.round(np.linspace(0.005, 0.05, 8), 4)),
    omega_r_values=list(np.round(np.linspace(0.03, 0.50, 8), 4)),
    chi_C_values=list(np.round(np.linspace(-0.9, 0.9, 8), 3)),
    omega_C_values=list(np.round(np.linspace(0.1, 0.9, 8), 3)),
)


def emit(name, y_axis, y_rng, x_axis, x_rng, disc):
    y_min, y_max, n_y = y_rng
    x_min, x_max, n_x = x_rng
    cfg = dict(BASE, y_axis=y_axis, y_min=y_min, y_max=y_max, n_y=n_y,
               x_axis=x_axis, x_min=x_min, x_max=x_max, n_x=n_x)
    (p1, v1), (p2, v2) = disc
    assert {y_axis, x_axis, p1, p2} == set(KNOBS), "axes + disc must cover all 4 knobs"
    sweeps, idx = [], 0
    for i1, a in enumerate(v1):
        for i2, b in enumerate(v2):
            sweeps.append({"index": idx, p1: float(a), p2: float(b),
                           f"i_{p1}": i1, f"i_{p2}": i2})
            idx += 1
    json.dump(dict(config=cfg, sweeps=sweeps), open(name, "w"), indent=2)
    print(f"{name}: {len(sweeps)} sweeps, axes ({y_axis} x {x_axis}) "
          f"{n_y}x{n_x}, discrete {p1} x {p2}")


emit("grid_chiC.json",
     "omega_C", (0.02, 0.90, 100), "chi_C", (-0.95, 0.95, 100),
     [("omega_r", BASE["omega_r_values"]), ("gamma", BASE["gamma_values"])])

emit("grid_omegar.json",
     "omega_C", (0.02, 0.90, 100), "omega_r", (0.02, 0.55, 100),
     [("chi_C", BASE["chi_C_values"]), ("gamma", BASE["gamma_values"])])

emit("grid_gammaomegar.json",
     "gamma", (0.005, 0.05, 100), "omega_r", (0.02, 0.55, 100),
     [("omega_C", BASE["omega_C_values"]), ("chi_C", BASE["chi_C_values"])])
