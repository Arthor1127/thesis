from collections import OrderedDict

import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.qubits.circle_transmon_squid import CircTransmonSQUID

# junction inductances (nH) - swept per EPR_C_EXTRACTION.md
L_TRANSMON_NH = 10.0   # single JJ to ground (doc's L_T) - jj_options.L_j
L_PROBE_NH = 15.0      # SQUID branch a (doc's L_probe) - squid_options.jj_a_options.L_j

# SQUID branch b is electrically bridged with continuous metal (export_mask=True
# below), not a real junction, per EPR_C_EXTRACTION.md Sec 4.1: one arm carries
# the probe inductor, the other is shorted so the loop is a single-inductor
# circuit. L_j/C_j here are unused electrically - qiskit_metal's option schema
# just needs numeric values for the (purely annotative) junction row.
L_BRIDGED_JJ_B_NH = 8.0
CJ_BRIDGED_JJ_B_FF = 1.8

# junction capacitances (fF) - fixed, not swept
CJ_TRANSMON_FF = 2.0
CJ_PROBE_FF = 1.5

# SQDMetal's create_port_JosephsonJunction addresses junctions on a component
# positionally (by insertion order into design.qgeometry.tables['junction']),
# not by name. CircTransmonSQUID emits: single JJ, then SQUID branch a, then b.
QUBIT_COMPONENT_NAME = 'Q1_SQUID'
JJ_INDEX_SINGLE = 0     # transmon junction -> L_TRANSMON_NH port
JJ_INDEX_SQUID_A = 1    # SQUID branch a -> L_PROBE_NH port
JJ_INDEX_SQUID_B = 2    # SQUID branch b -> bridged with metal, no port created

# resonator parameters
RES_LENGTH_MM = 3.74575

# chip geometry: size_x/size_y/center_x/center_y are auto-computed from the
# drawn geometry's bounding box (see _design_bounds_mm / CHIP_MARGIN_FACTOR
# below) unless explicitly overridden. size_z is fixed.
CHIP_SIZE_Z_UM = 500
CHIP_MARGIN_FACTOR = 4.0   # chip footprint = design bounding box * this factor


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


def build_design(L_transmon_nH=L_TRANSMON_NH, L_probe_nH=L_PROBE_NH,
                  res_length_mm=RES_LENGTH_MM,
                  chip_size_x_mm=None, chip_size_y_mm=None,
                  chip_center_x_mm=None, chip_center_y_mm=None):
    """
    Build the circular transmon + SQUID design coupled to a meandered
    readout resonator running between an open and a short termination.

    The SQUID's branch b is bridged with continuous metal (export_mask=True)
    so the loop is a single-inductor circuit with branch a as the sole
    junction - see EPR_C_EXTRACTION.md Sec 4.1.

    The chip is sized AFTER the geometry is drawn, from the geometry's own
    bounding box * CHIP_MARGIN_FACTOR, centered on that bounding box - so the
    chip grows/shrinks around wherever the components already are instead of
    the components being placed relative to a fixed chip size. Any of the
    chip_* kwargs can be passed to override the corresponding auto-computed
    value.
    """
    design = designs.DesignPlanar({}, overwrite_enabled=True)

    design.variables['cpw_width'] = '20um'
    design.variables['cpw_gap'] = '12.25um'

    # --- Qubit (circular transmon with SQUID) ---

    qubit_options = Dict(
        pos_x='0um', pos_y='0um', orientation='0', chip='main',
        cpw_width='cpw_width', cpw_gap='cpw_gap',

        pad_radius='90um',
        pad_gap='10um',

        jj_options=Dict(
            jj_width='5um',
            jj_angle='0',
            L_j=f'{L_transmon_nH}nH',
            C_j=f'{CJ_TRANSMON_FF}fF',
            export_mask=False,
            jj_sim_gap='3um',
        ),

        squid_options=Dict(
            theta_2='180',
            delta_angle='7.5',
            g1='4um',

            l1='600um',
            w1='6um',
            l2='600um',

            jj_a_options=Dict(
                jj_width='4um', L_j=f'{L_probe_nH}nH', C_j=f'{CJ_PROBE_FF}fF',
                export_mask=False, jj_sim_gap='3um',
            ),
            jj_b_options=Dict(
                jj_width='6um', L_j=f'{L_BRIDGED_JJ_B_NH}nH', C_j=f'{CJ_BRIDGED_JJ_B_FF}fF',
                export_mask=True, jj_sim_gap='3um',
            ),
        ),
    )
    qubit_1 = CircTransmonSQUID(design, QUBIT_COMPONENT_NAME, options=qubit_options)

    # --- Readout resonator terminations ---

    launch_point1 = OpenToGround(design, 'launch_point1', options=dict(
        pos_x='-275um', pos_y='-850um', orientation='0',
        termination_gap='cpw_gap', gap='cpw_gap', width='cpw_width'))

    launch_point2 = ShortToGround(design, 'launch_point2', options=dict(
        pos_x='-150um', pos_y='50um', orientation='0',
        termination_gap='cpw_gap', gap='cpw_gap', width='cpw_width'))

    jogs = OrderedDict()
    jogs[0] = ["L", "100um"]
    jogs[1] = ["L", "700um"]

    readout_1 = RouteMeander(design, 'readout_1', options=dict(
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
    print(f"L_TRANSMON_NH={L_TRANSMON_NH}, L_PROBE_NH={L_PROBE_NH}, RES_LENGTH_MM={RES_LENGTH_MM}")
    d = build_design()
    print("Design built OK. Components:", list(d.components.keys()))
