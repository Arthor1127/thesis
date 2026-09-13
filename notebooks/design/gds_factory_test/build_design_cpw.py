"""Bare meandered CPW readout resonator, stripped from
``squid_array_coupling/build_design.py`` with the SquidArray removed -- just
the OpenToGround / ShortToGround / RouteMeander chain and the auto-sized
chip. Used to test GroundPlaneCheese against a genuinely non-rectangular
keepout shape (a fillet-jogged meander), not just a component's own straight
bounding box.
"""

from collections import OrderedDict

import numpy as np
from qiskit_metal import designs
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround

# ---- readout resonator ----
RES_LENGTH_MM = 3.74575

# ---- chip geometry ----
# chip footprint = design bounding box * this factor (see squid_array_coupling
# for why 2.5, not a larger margin).
CHIP_SIZE_Z_UM = 500
CHIP_MARGIN_FACTOR = 2.5


def _design_bounds_mm(design):
    """Bounding box (minx, miny, maxx, maxy), in mm, across every drawn
    qgeometry table. Must be called AFTER all components are built."""
    all_bounds = []
    for _, table in design.qgeometry.tables.items():
        if len(table) == 0:
            continue
        all_bounds.append(table.total_bounds)
    all_bounds = np.array(all_bounds)
    return (all_bounds[:, 0].min(), all_bounds[:, 1].min(),
            all_bounds[:, 2].max(), all_bounds[:, 3].max())


def build_design(res_length_mm=RES_LENGTH_MM,
                 chip_size_x_mm=None, chip_size_y_mm=None,
                 chip_center_x_mm=None, chip_center_y_mm=None):
    """Build the meandered-resonator-only design (no SQUID array)."""
    design = designs.DesignPlanar({}, overwrite_enabled=True)

    design.variables['cpw_width'] = '20um'
    design.variables['cpw_gap'] = '12.25um'

    OpenToGround(design, 'launch_point1', options=dict(
        pos_x='-275um', pos_y='-850um', orientation='0',
        termination_gap='cpw_gap', gap='cpw_gap', width='cpw_width'))

    ShortToGround(design, 'launch_point2', options=dict(
        pos_x='-150um', pos_y='50um', orientation='0',
        termination_gap='cpw_gap', gap='cpw_gap', width='cpw_width'))

    jogs = OrderedDict()
    jogs[0] = ["L", "100um"]
    jogs[1] = ["L", "700um"]

    RouteMeander(design, 'readout_1', options=dict(
        pin_inputs=dict(
            start_pin=dict(component='launch_point1', pin='open'),
            end_pin=dict(component='launch_point2', pin='short'),
        ),
        fillet='40um',
        lead=dict(
            start_straight='10um',
            end_straight='1180um',
            end_jogged_extension=jogs,
        ),
        total_length=f'{res_length_mm}mm',
    ))

    # --- Size the chip from the drawn geometry's own bounding box ---
    minx, miny, maxx, maxy = _design_bounds_mm(design)
    bounds_center_x, bounds_center_y = (minx + maxx) / 2, (miny + maxy) / 2
    bounds_size_x, bounds_size_y = (maxx - minx), (maxy - miny)

    chip_size_x_mm = bounds_size_x * CHIP_MARGIN_FACTOR if chip_size_x_mm is None else chip_size_x_mm
    chip_size_y_mm = bounds_size_y * CHIP_MARGIN_FACTOR if chip_size_y_mm is None else chip_size_y_mm
    chip_center_x_mm = bounds_center_x if chip_center_x_mm is None else chip_center_x_mm
    chip_center_y_mm = bounds_center_y if chip_center_y_mm is None else chip_center_y_mm

    design.chips.main.size.size_x = f'{chip_size_x_mm}mm'
    design.chips.main.size.size_y = f'{chip_size_y_mm}mm'
    design.chips.main.size.size_z = f'{CHIP_SIZE_Z_UM}um'
    design.chips.main.size.center_x = f'{chip_center_x_mm}mm'
    design.chips.main.size.center_y = f'{chip_center_y_mm}mm'

    return design


if __name__ == "__main__":
    d = build_design()
    mnx, mny, mxx, mxy = _design_bounds_mm(d)
    print("Design built OK. Components:", list(d.components.keys()))
    print(f"Design bounds (mm): X=[{mnx:.4f},{mxx:.4f}]  Y=[{mny:.4f},{mxy:.4f}]")
    print(f"Chip size (mm): {d.chips.main.size.size_x} x {d.chips.main.size.size_y}")
