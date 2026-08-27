import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.lumped.cap_n_interdigital import CapNInterdigital

# Capacitor geometry (interdigital finger cap) - swept per lumped.ipynb
CAP_FINGER_COUNT = 5
CAP_GAP_UM = 6.0          # gap between fingers - the main capacitance knob

# Capacitor geometry - fixed
CAP_WIDTH_UM = 10.0
CAP_FINGER_LENGTH_UM = 20.0
CAP_DISTANCE_UM = 50.0
CAP_NS_WIDTH_UM = 10.0
CAP_NS_GAP_UM = 6.0
CAP_GAP_GROUND_UM = 6.0

# Meander inductor - total_length is the main inductance knob
RES_TOTAL_LENGTH_MM = 5.0

# Meander geometry - fixed. spacing/fillet/lead_straight must satisfy
# fillet < ~spacing/2 and fillet < ~lead_straight/2, or corners near/at the
# capacitor and inside the meander itself will stay square (qiskit_metal
# logs a "short segments" warning per corner - see lumped.ipynb for a worked
# example of the zero-meander and unfilleted-corner failure modes this
# avoids).
MEANDER_FILLET_UM = 39.9
MEANDER_LEAD_STRAIGHT_UM = 200.0
MEANDER_SPACING_UM = 80.0
MEANDER_ASYMMETRY_UM = 850.0

# chip geometry: size_x/size_y/center_x/center_y are auto-computed from the
# drawn geometry's bounding box (see _design_bounds_mm / CHIP_MARGIN_FACTOR
# below) unless explicitly overridden. size_z is fixed.
CHIP_SIZE_Z_UM = 500
CHIP_MARGIN_FACTOR = 3.0   # chip footprint = design bounding box * this factor


def _design_bounds_mm(design):
    """Bounding box (minx, miny, maxx, maxy), in mm, across every drawn
    qgeometry table. Must be called AFTER all components are built."""
    all_bounds = []
    for table_name, table in design.qgeometry.tables.items():
        if len(table) == 0:
            continue
        all_bounds.append(table.total_bounds)  # [minx, miny, maxx, maxy], mm
    all_bounds = np.array(all_bounds)
    return (all_bounds[:, 0].min(), all_bounds[:, 1].min(),
            all_bounds[:, 2].max(), all_bounds[:, 3].max())


def build_design(finger_count=CAP_FINGER_COUNT, cap_gap_um=CAP_GAP_UM,
                  res_length_mm=RES_TOTAL_LENGTH_MM,
                  chip_size_x_mm=None, chip_size_y_mm=None,
                  chip_center_x_mm=None, chip_center_y_mm=None):
    """
    Build the lumped-element LC resonator from lumped.ipynb: an interdigital
    finger capacitor (`CapNInterdigital`) whose two pins (`north_end`,
    `south_end`) are closed into a loop by a meandered inductor
    (`RouteMeander`). The structure floats on its own - no feedline, no
    ground short, no junction - so eigenmode extraction gives the bare LC
    resonance directly.

    The chip is sized AFTER the geometry is drawn, from the geometry's own
    bounding box * CHIP_MARGIN_FACTOR, centered on that bounding box - so the
    chip grows/shrinks around wherever the components already are instead of
    the components being placed relative to a fixed chip size. Any of the
    chip_* kwargs can be passed to override the corresponding auto-computed
    value.
    """
    design = designs.DesignPlanar({}, overwrite_enabled=True)

    # --- Interdigital capacitor ---

    cap_options = Dict(
        finger_count=str(int(finger_count)),
        cap_gap=f'{cap_gap_um}um',
        cap_width=f'{CAP_WIDTH_UM}um',
        finger_length=f'{CAP_FINGER_LENGTH_UM}um',
        cap_distance=f'{CAP_DISTANCE_UM}um',
        north_width=f'{CAP_NS_WIDTH_UM}um',
        north_gap=f'{CAP_NS_GAP_UM}um',
        south_width=f'{CAP_NS_WIDTH_UM}um',
        south_gap=f'{CAP_NS_GAP_UM}um',
        cap_gap_ground=f'{CAP_GAP_GROUND_UM}um',
    )
    cap = CapNInterdigital(design, 'cap', options=cap_options)

    # --- Meandered inductor closing the loop across the capacitor ---

    inductor = RouteMeander(design, 'inductor', options=Dict(
        total_length=f'{res_length_mm}mm',
        trace_width='cpw_width',
        fillet=f'{MEANDER_FILLET_UM}um',
        pin_inputs=Dict(
            start_pin=Dict(component='cap', pin='north_end'),
            end_pin=Dict(component='cap', pin='south_end'),
        ),
        lead=Dict(
            start_straight=f'{MEANDER_LEAD_STRAIGHT_UM}um',
            end_straight=f'{MEANDER_LEAD_STRAIGHT_UM}um',
        ),
        meander=Dict(
            spacing=f'{MEANDER_SPACING_UM}um',
            asymmetry=f'{MEANDER_ASYMMETRY_UM}um',
        ),
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
    print(f"CAP_FINGER_COUNT={CAP_FINGER_COUNT}, CAP_GAP_UM={CAP_GAP_UM}, "
          f"RES_TOTAL_LENGTH_MM={RES_TOTAL_LENGTH_MM}")
    d = build_design()
    print("Design built OK. Components:", list(d.components.keys()))
