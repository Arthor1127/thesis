"""Geometry-only check: build the design and render it, with NO meshing
and NO Palace involvement, so the geometry itself can be reviewed cheaply
before committing to a full simulation."""
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import build_design, _design_bounds_mm

design = build_design()
design.rebuild()

minx, miny, maxx, maxy = _design_bounds_mm(design)
print(f"Design bounds (mm): X=[{minx:.4f},{maxx:.4f}]  Y=[{miny:.4f},{maxy:.4f}]")
print(f"Chip size (mm): {design.chips.main.size.size_x} x {design.chips.main.size.size_y}, "
      f"centered at ({design.chips.main.size.center_x}, {design.chips.main.size.center_y})")
print(f"Components: {list(design.components.keys())}")

from qiskit_metal import MetalGUI
gui = MetalGUI(design)
gui.rebuild()
gui.autoscale()
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "geometry_test.png")
gui.screenshot(display=False, name=out_path)
print(f"Saved geometry screenshot to {out_path}")
