"""Write grid.json for the stanene (B_x, epsilon_site) gap map.

Matches Analisis_gaps_cris.ipynb: mesh x mesh points, B in [B_min, B_max],
E in [E_min, E_max], one array task per B row.
"""
import json

GRID = dict(
    mesh=50,
    B_min=0.0, B_max=1.0,
    E_min=0.0, E_max=0.6,
    gap2d_mallado=240,
    gap1d_mallado=240,
    ribbon_num_cells=35,
)

if __name__ == "__main__":
    json.dump(GRID, open("grid.json", "w"), indent=2)
    m = GRID["mesh"]
    print("grid.json:", GRID, "->", m, "row tasks,", m * m, "points total")
