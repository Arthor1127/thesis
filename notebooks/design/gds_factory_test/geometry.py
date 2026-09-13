"""
Minimal Qiskit Metal test design: one transmon qubit with a meandered
readout resonator terminated in an open. Exports to GDS with Qiskit
Metal's own cheesing disabled, ready to feed straight into
vortex_trap_fill_workflow.py.

No GUI/Qt needed -- runs headless. If you want to look at it first,
use `qm.view(design)` (see comment at the bottom).
"""

from qiskit_metal import designs
from qiskit_metal.qlibrary.qubits.transmon_pocket import TransmonPocket
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander

design = designs.DesignPlanar()
design.overwrite_enabled = True
design.chips.main.size.size_x = "5mm"
design.chips.main.size.size_y = "5mm"

q1 = TransmonPocket(
    design,
    "Q1",
    options=dict(
        pad_width="425um",
        pocket_height="650um",
        connection_pads=dict(readout=dict(loc_W=+1, loc_H=+1, pad_width="200um")),
    ),
)

otg = OpenToGround(
    design,
    "open_to_ground",
    options=dict(pos_x="1.75mm", pos_y="0um", orientation="0"),
)

readout = RouteMeander(
    design,
    "readout",
    options=dict(
        total_length="6mm",
        hfss_wire_bonds=True,
        fillet="90um",
        lead=dict(start_straight="100um"),
        pin_inputs=dict(
            start_pin=dict(component="Q1", pin="readout"),
            end_pin=dict(component="open_to_ground", pin="open"),
        ),
    ),
)

# Disable Qiskit Metal's own cheesing -- the qpdk/gdsfactory pipeline
# does the fill instead, so we don't want two overlapping hole patterns.
design.renderers.gds.options.cheese["view_in_file"]["main"] = {1: False}
design.renderers.gds.options.no_cheese["view_in_file"]["main"] = {1: False}

design.renderers.gds.export_to_gds("my_design.gds")
print("Wrote my_design.gds")
print(
    "Main-chip metal is on GDS layer 1, datatype 0 by default -- "
    "matches METAL_LAYER = (1, 0) in vortex_trap_fill_workflow.py."
)

# To preview without a GUI:
# import qiskit_metal as qm
# fig = qm.view(design)
# fig.savefig("my_design_preview.png")