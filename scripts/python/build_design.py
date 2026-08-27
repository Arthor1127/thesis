"""
Builds the schizo-resonator qiskit-metal design.

Topology: open-open (otg1, otg2 both OpenToGround) half-wave resonator,
with an optional single grounding strap at a swept position along its length.

--- Length calibration (5 GHz fundamental) ---
Reference data point (measured, quarter-wave open-short topology):
    f1_ref = 7.980849 GHz at total_length = 3.7mm  =>  f1 = v/(4L)
    v = 4 * 3.7mm * 7.980849 GHz = 1.1812e8 m/s

Our design is open-open (half-wave): f1 = v/(2L). Assuming the SAME phase
velocity (same substrate/cross-section/trace geometry), solving for L at
f1 = 5 GHz:
    L = v / (2 * 5 GHz) = 11.812 mm

CAVEAT: this rescales across two different topologies (quarter-wave -> half-
wave), assuming phase velocity is topology-independent (a reasonable but
unverified assumption). The ground_pos_mm=None (ungrounded) baseline run in
the sweep exists specifically to measure the TRUE f1 of this exact geometry
so TOTAL_LENGTH_MM can be corrected if the baseline comes back off-target.
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["PMIX_MCA_gds"] = "hash"

from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.meandered_grounded import RouteMeanderGrounded
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround

# Derived length for a 5 GHz half-wave fundamental - see docstring above.
TOTAL_LENGTH_MM = 14.9669
GROUND_WIDTH_UM = 10
# Keep the strap this far from either physical end of the meander path,
# since a strap placed exactly at an endpoint is not physically meaningful.
EDGE_MARGIN_MM = 0.5


def build_design(ground_pos_mm=None, total_length_mm=TOTAL_LENGTH_MM):
    """
    Build the design.

    ground_pos_mm : float or None
        Position (mm) along the resonator path for a single grounding
        strap. If None, builds the UNGROUNDED baseline (plain RouteMeander,
        no strap) - this is the reference point used to validate/calibrate
        TOTAL_LENGTH_MM against the true measured f1.
    total_length_mm : float
        Physical path length of the meander, in mm.
    """
    design = designs.DesignPlanar({}, overwrite_enabled=True)

    design.chips.main.size.size_x = '4.8mm'
    design.chips.main.size.size_y = '2.4mm'
    design.chips.main.size.size_z = '500um'
    design.chips.main.size.center_x = '0mm'
    design.chips.main.size.center_y = '-1mm'

    design.variables['cpw_width'] = '10 um'
    design.variables['cpw_gap'] = '6 um'

    x1, y1 = '-2mm', '0mm'
    launch_options1 = dict(chip='main', pos_x=x1, pos_y=y1, orientation='360',
                            lead_length='30um', pad_height='103um',
                            pad_width='103um', pad_gap='60um')
    LP1 = LaunchpadWirebond(design, 'LP1', options=launch_options1)

    x2 = '2mm'
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

    # Both ends OPEN (this is what makes it a half-wave resonator, and
    # what the S21 driven sweep request means by "both ends open")
    otg1 = OpenToGround(design, 'otg1', options=dict(
        chip='main', pos_x='-0.2mm', pos_y='-40um', orientation=180))
    otg2 = OpenToGround(design, 'otg2', options=dict(
        chip='main', pos_x='0mm', pos_y='-1.35mm', orientation=-90))

    common_kwargs = dict(
        trace_width='10um',
        trace_gap='6um',
        total_length=f'{total_length_mm}mm',
        hfss_wire_bonds=False,
        fillet='99.9 um',
        lead=dict(start_straight='300um'),
        pin_inputs=Dict(
            start_pin=Dict(component='otg1', pin='open'),
            end_pin=Dict(component='otg2', pin='open')),
    )

    if ground_pos_mm is None:
        # Baseline: no grounding strap at all - plain half-wave resonator
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
                width=f'{GROUND_WIDTH_UM}um'),
        ))

    return design


def sweep_positions_mm(n_points=9, total_length_mm=TOTAL_LENGTH_MM,
                        edge_margin_mm=EDGE_MARGIN_MM):
    """Evenly spaced grounding-strap positions spanning the full resonator
    length (minus edge margins). Does NOT include the None (baseline) case
    - add that separately in the caller."""
    import numpy as np
    return [round(float(x), 4) for x in
            np.linspace(edge_margin_mm, total_length_mm - edge_margin_mm, n_points)]


if __name__ == "__main__":
    print(f"TOTAL_LENGTH_MM = {TOTAL_LENGTH_MM}")
    print(f"Sweep positions (mm): {sweep_positions_mm()}")

    d_baseline = build_design(ground_pos_mm=None)
    print("Baseline (ungrounded) design built OK. Components:", list(d_baseline.components.keys()))

    d_grounded = build_design(ground_pos_mm=5.9)
    print("Grounded design built OK. Components:", list(d_grounded.components.keys()))
