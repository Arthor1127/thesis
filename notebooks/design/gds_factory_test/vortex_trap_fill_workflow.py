"""
Vortex-trap hole fill for a Qiskit Metal chip, applied via gdsfactory/qpdk.

Places holes on a triangular (equilateral-triangle) lattice rather than a
plain square grid, by choosing gap/stagger values that turn qpdk's general
oblique-lattice fill into a regular triangular one.

WORKFLOW
--------
1. Finish your design in Qiskit Metal as normal.
2. Turn OFF Qiskit Metal's own "cheesing" before export, so you don't end
   up with two overlapping hole patterns:

       design.renderers.gds.options.cheese["view_in_file"]["main"] = {1: False}
       design.renderers.gds.export_to_gds("my_design.gds")

   (If your Qiskit Metal version still writes a cheese/no_cheese layer
   regardless, that's fine -- just note its layer/datatype and strip it
   out below with `chip.remove_layers([...])` before filling.)

3. Run this script (adjust the CONFIG section below to match your design).
4. Fabricate / inspect `my_design_filled.gds`.

Notes:
- gdsfactory's `import_gds` preserves the exact (layer, datatype) numbers
  from the source GDS -- it does NOT remap them to qpdk's own layer names.
  Open the exported GDS in KLayout (or check the printed layer list below)
  to find the actual layer/datatype numbers for your ground plane / metal
  traces before setting METAL_LAYER.
- Vortex-trap fill is a fabrication-mask step, not something you'd feed
  back into Qiskit Metal for EM simulation -- keep the unfilled and filled
  GDS files separate.
"""

import math
import gdsfactory as gf
from qpdk import PDK
from qpdk.utils import fill_magnetic_vortices


PDK.activate()

INPUT_GDS = "my_design.gds"  # exported from Qiskit Metal
OUTPUT_GDS = "my_design_filled.gds"
METAL_LAYER = (1, 0)  # <-- your ground plane / metal layer, datatype
FILL_LAYER = (66, 0)  # <-- free layer/datatype for the holes
KEEPOUT_UM = 50.0  # clearance to keep free of holes around circuit metal
PITCH_UM = 30.0  # center-to-center hole spacing (tune to your film/design)
HOLE_SIZE_UM = (8.0, 8.0)  # (width, height) of each hole

chip = gf.import_gds(INPUT_GDS)
print("Layers found:", sorted(chip.layers))

# Triangular-lattice gap/stagger for the given pitch (row spacing = triangle height)
gap = (PITCH_UM - HOLE_SIZE_UM[0], PITCH_UM * math.sqrt(3) / 2 - HOLE_SIZE_UM[1])
stagger = (-PITCH_UM / 2, 0.0)
assert min(gap) > 0, "hole_size too large for the chosen pitch"

filled_chip = fill_magnetic_vortices(
    component=chip,
    rectangle_size=HOLE_SIZE_UM,
    gap=gap,
    stagger=stagger,
    exclude_layers=[(METAL_LAYER, KEEPOUT_UM)],
    fill_layer=FILL_LAYER,
)
filled_chip.write_gds(OUTPUT_GDS)
print(f"Wrote {OUTPUT_GDS}")

# Optional: inspect interactively
# filled_chip.show()
