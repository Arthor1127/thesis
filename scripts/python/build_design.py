from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround

TOTAL_LENGTH_MM = 7.48345

def build_design(total_length_mm=TOTAL_LENGTH_MM):

    total_length_mm = TOTAL_LENGTH_MM
    design = designs.DesignPlanar({}, overwrite_enabled=True)
    design.chips.main.size.size_x = '4.0mm'
    design.chips.main.size.size_y = '3mm'
    design.chips.main.size.size_z = '500um' # tamaño del sustrato y de la capa de aire a simular
    design.chips.main.size.center_x = '0mm'
    design.chips.main.size.center_y = '-0.5mm'

    design.variables['cpw_width'] = '10 um'
    design.variables['cpw_gap'] = '6 um'

    x1, y1 = '-1mm', '0mm'
    launch_options1 = dict(chip='main', pos_x=x1, pos_y=y1, orientation='360',
                        lead_length='30um', pad_height='103um',
                        pad_width='103um', pad_gap='60um')
    LP1 = LaunchpadWirebond(design, 'LP1', options=launch_options1)

    x2 = '1mm'
    launch_options2 = dict(chip='main', pos_x=x2, pos_y=y1, orientation='180',
                            lead_length='30um', pad_height='103um',
                            pad_width='103um', pad_gap='60um')
    LP2 = LaunchpadWirebond(design, 'LP2', options=launch_options2)

    TL = RoutePathfinder(design, 'TL', options=dict(
        chip='main', trace_width='10um', trace_gap='6um', fillet='90um',
        hfss_wire_bonds=True, lead=dict(end_straight='0.1mm'),
        pin_inputs=Dict(
            start_pin=Dict(component='LP1', pin='tie'),
            end_pin=Dict(component='LP2', pin='tie'))))
    otg1 = OpenToGround(design, 'otg1', options=dict(
        chip='main', pos_x='-0.2mm', pos_y='-40um', orientation=180))
    stg2 = ShortToGround(design, 'stg2', options=dict(
        chip='main', pos_x='0mm', pos_y='-1.25mm', orientation=0))

    common_kwargs = dict(
        trace_width='10um',
        trace_gap='6um',
        total_length=f'{total_length_mm}mm',
        hfss_wire_bonds=False,
        fillet='99.9 um',
        lead=dict(start_straight='10um'),
        meander = Dict(spacing = "200um", assymetry = "1mm"),
        pin_inputs=Dict(
            start_pin=Dict(component='otg1', pin='open'),
            end_pin=Dict(component='stg2', pin='short')),
    )

    res1 = RouteMeander(design, 'resonator1', Dict(**common_kwargs))

    return design