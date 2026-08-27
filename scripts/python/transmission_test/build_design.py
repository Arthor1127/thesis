from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.meandered_grounded import RouteMeanderGrounded
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight
from qiskit_metal import MetalGUI

# resonator parameters
TOTAL_LENGTH_MM = 5.92905
GROUND_WIDTH_UM = 10
GROUND_OVERLAP_UM = 4
EDGE_MARGIN_MM = 0.5

# chip geometry (make it as tight as possible)
CHIP_SIZE_X_MM = 3.0
CHIP_SIZE_Y_MM = 2.0
CHIP_SIZE_Z_UM = 500
CHIP_CENTER_X_MM = 0.0
CHIP_CENTER_Y_MM = -0.381

def build_design(ground_pos_mm=None, total_length_mm=TOTAL_LENGTH_MM,
                  chip_size_x_mm=None, chip_size_y_mm=None,
                  chip_center_x_mm=None, chip_center_y_mm=None,
                  ground_overlap_um=None):
    """
    Build the quarter-wave resonator design WITH a feedline for driven simulation.
    """
    chip_size_x_mm = CHIP_SIZE_X_MM if chip_size_x_mm is None else chip_size_x_mm
    chip_size_y_mm = CHIP_SIZE_Y_MM if chip_size_y_mm is None else chip_size_y_mm
    chip_center_x_mm = CHIP_CENTER_X_MM if chip_center_x_mm is None else chip_center_x_mm
    chip_center_y_mm = CHIP_CENTER_Y_MM if chip_center_y_mm is None else chip_center_y_mm
    ground_overlap_um = GROUND_OVERLAP_UM if ground_overlap_um is None else ground_overlap_um
    
    design = designs.DesignPlanar({}, overwrite_enabled=True)

    design.chips.main.size.size_x = f'{chip_size_x_mm}mm'
    design.chips.main.size.size_y = f'{chip_size_y_mm}mm'
    design.chips.main.size.size_z = f'{CHIP_SIZE_Z_UM}um'
    design.chips.main.size.center_x = f'{chip_center_x_mm}mm'
    design.chips.main.size.center_y = f'{chip_center_y_mm}mm'

    design.variables['cpw_width'] = '20 um'
    design.variables['cpw_gap'] = '12.25 um'

    fillet='99.99um'

    cpw_options = Dict(
        lead=Dict(
            start_straight='100um',
            end_straight='100um'),
        fillet=fillet
        )

    main_branch_height_um= 50
    main_branch_half_width_um = 800
    
    pad_1_options = dict(
            pos_x = f"-{main_branch_half_width_um}um",
            pos_y = f"{main_branch_height_um}um",
            orientation = "0"
        )

    pad_2_options = dict(
            pos_x = f"{main_branch_half_width_um}um",
            pos_y = f"{main_branch_height_um}um",
            orientation = "180"
        )
    
    # Renamed to LP1
    pad_1 = LaunchpadWirebond(
        design, 
        "LP1",
        options = pad_1_options
    )
    # Renamed to LP2
    pad_2 = LaunchpadWirebond(
        design, 
        "LP2",
        options = pad_2_options
    )
    main_branch_options = {
        "pin_inputs": {
            "start_pin": {
                "component": "LP1", 
                "pin": "tie"
                },
            "end_pin": {
                "component": "LP2", 
                "pin": "tie"
                }
        }
    }

    # Renamed to TL
    main_branch = RouteStraight(design, "TL", main_branch_options)
    # ---------------------------------------------------------------------

    # --- Resonator Components ---
    
    otg1 = OpenToGround(design, 'otg1', options=dict(
        chip='main', pos_x='-0.5mm', pos_y='-10um', orientation=180,
        width='20um', gap='12.25um', termination_gap='12.25um'))
        
    sg1 = ShortToGround(design, 'sg1', options=dict(
        chip='main', pos_x='0mm', pos_y='-0.91mm', orientation=-90,
        width='20um', gap='12.25um'))

    common_kwargs = dict(
        trace_width='20um',
        trace_gap='12.25um',
        total_length=f'{total_length_mm}mm',
        hfss_wire_bonds=True,
        fillet='99.9 um',
        lead=dict(start_straight='100um', end_straight = "20um"),
        pin_inputs=Dict(
            start_pin=Dict(component='sg1', pin='short'),
            end_pin=Dict(component='otg1', pin='open')),
    )

    if ground_pos_mm is None:
        res1 = RouteMeander(design, 'resonator1', Dict(**common_kwargs))
    else:
        if not (EDGE_MARGIN_MM <= ground_pos_mm <= total_length_mm - EDGE_MARGIN_MM):
            raise ValueError(
                f"ground_pos_mm={ground_pos_mm} out of valid range "
                f"[{EDGE_MARGIN_MM}, {total_length_mm - EDGE_MARGIN_MM}]"
            )
        res1 = RouteMeanderGrounded(design, 'resonator1', Dict(
            **common_kwargs,
            ground_straps=dict(
                positions=[f'{ground_pos_mm}mm'],
                width=f'{GROUND_WIDTH_UM}um',
                overlap=f'{ground_overlap_um}um'),
        ))

    return design

if __name__ == "__main__":
    print(f"TOTAL_LENGTH_MM = {TOTAL_LENGTH_MM}")
    d_grounded = build_design(ground_pos_mm=2.0)
    print("Grounded design built OK. Components:", list(d_grounded.components.keys()))