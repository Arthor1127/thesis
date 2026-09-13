"""Exports the bare shorted quarter-wave CPW resonator geometry used in
run_bareres_eigenmode_isolation.py (parametric_coupling_chip/) to GDS, for
the ANSYS cross-check requested in this directory's HANDOFF.md.

Geometry-only (qiskit_metal, no PALACE) -- safe to run locally, matching this
repo's standing rule that only actual FEM/MPI solves go to the cluster.

Run: conda activate qcg-quantum-design && python export_geometry.py
Writes: geometry/bare_resonator.gds
"""
import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder

CPW_WIDTH_UM = 20.0
CPW_GAP_UM = 11.5
TARGET_RES_LENGTH_MM = 3.77259   # same as run_bareres_eigenmode_isolation.py / readout_resonators.py
half_span_um = TARGET_RES_LENGTH_MM * 1e3 / 2.0
chip_size_x_mm = (TARGET_RES_LENGTH_MM * 1e3 + 800.0) / 1e3

design = designs.DesignPlanar({}, overwrite_enabled=True)
design.chips.main.size.size_x = f'{chip_size_x_mm}mm'
design.chips.main.size.size_y = '2mm'
design.chips.main.size.size_z = '500um'
design.chips.main.size.center_x = '0mm'
design.chips.main.size.center_y = '0mm'
design.variables['cpw_width'] = f'{CPW_WIDTH_UM}um'
design.variables['cpw_gap'] = f'{CPW_GAP_UM}um'

x1 = f'-{half_span_um}um'
x2 = f'{half_span_um}um'
OpenToGround(design, 'OTG1', options=dict(
    chip='main', pos_x=x1, pos_y='0um', orientation='180', width='cpw_width'))
ShortToGround(design, 'STG', options=dict(
    chip='main', pos_x=x2, pos_y='0um', orientation='0', width='cpw_width'))
RoutePathfinder(design, 'RES', options=dict(
    chip='main', trace_width='cpw_width', trace_gap='cpw_gap', fillet='50um',
    pin_inputs=Dict(start_pin=Dict(component='OTG1', pin='open'),
                    end_pin=Dict(component='STG', pin='short'))))
design.rebuild()

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'geometry', 'bare_resonator.gds')
design.renderers.gds.export_to_gds(out_path)
print(f"wrote {out_path}")
print(f"chip size: {design.chips.main.size.size_x} x {design.chips.main.size.size_y}")
print(f"OTG1 (open end, where the lossy port goes) at x={x1}")
print(f"STG (real PEC short) at x={x2}")
