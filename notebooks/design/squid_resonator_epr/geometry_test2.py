"""Geometry-only check, take 2: plot the actual QGeometry polygons directly
via matplotlib (no GUI screenshot ambiguity), zoomed on the qubit/SQUID
region specifically."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import build_design, _design_bounds_mm

design = build_design()
design.rebuild()

minx, miny, maxx, maxy = _design_bounds_mm(design)
print(f"Design bounds (mm): X=[{minx:.4f},{maxx:.4f}]  Y=[{miny:.4f},{maxy:.4f}]")

fig, ax = plt.subplots(figsize=(10, 8))
for table_name, table in design.qgeometry.tables.items():
    if len(table) == 0:
        continue
    for geom in table.geometry:
        if geom.geom_type == "Polygon":
            xs, ys = geom.exterior.xy
            ax.fill(xs, ys, alpha=0.6, edgecolor="black", linewidth=0.3)
        elif geom.geom_type == "LineString":
            xs, ys = geom.xy
            ax.plot(xs, ys, linewidth=1)
        elif geom.geom_type == "MultiPolygon":
            for poly in geom.geoms:
                xs, ys = poly.exterior.xy
                ax.fill(xs, ys, alpha=0.6, edgecolor="black", linewidth=0.3)

ax.set_aspect("equal")
ax.set_xlabel("x (mm)")
ax.set_ylabel("y (mm)")
ax.set_title("Full design")
fig.savefig("geometry_test_full.png", dpi=150)
print("Saved geometry_test_full.png")

# Zoomed on the qubit/SQUID region
ax.set_xlim(-1.4, 0.2)
ax.set_ylim(-0.3, 0.3)
ax.set_title("Zoomed on qubit/SQUID region")
fig.savefig("geometry_test_qubit_zoom.png", dpi=150)
print("Saved geometry_test_qubit_zoom.png")
