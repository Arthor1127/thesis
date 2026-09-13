"""Geometry for validating what SQDMetal/Palace's `create_port_2_conds`
actually realizes for `impedance_C`, independent of any resonator.

WHY: the two-port hanger isolation test (hanger_isolation_design.py) gave
Qc(from kappa)/Qc_target = 0.08, 8.86, 2.67 -- wildly inconsistent, and at
the strongest coupling the resonator mode couldn't even be found near its
expected frequency. Leading hypothesis: `create_port_2_conds`'s impedance_C
is realized through some geometry-dependent conversion, the same way we
already know (and never resolved) PALACE/SQDMetal converts a target port
*resistance* into a geometry-dependent sheet impedance rather than the
nominal value directly. This test isolates that question completely: no
resonator, no mode-finding, just a plain shunt capacitor on a matched
through-line, with an exact closed-form S21(f) prediction to check the
DRIVEN S-parameter simulation against directly, at any frequency.

CIRCUIT: a lumped admittance Y=j*omega*C shunted to ground at one point on
an otherwise uniform, matched (Z0 on both sides) transmission line has the
standard two-port network-theory result (e.g. Pozar):
    S21 = 1 / (1 + Y*Z0/2)         S11 = -(Y*Z0/2) / (1 + Y*Z0/2)
This is EXACT for an ideal lumped shunt element -- no resonator, no
approximation, valid at every frequency. Comparing PALACE's simulated S21(f)
against this formula, at the SAME Cc values used in the hanger isolation
test (2.88/9.11/28.8 fF), directly reveals whether create_port_2_conds
realizes the requested capacitance -- and if not, by what factor.

GEOMETRY: identical through-line to hanger_isolation_design.py (OTG1--TL1--
[LineTee]--TL2--OTG2), but the LineTee's `second_end` pin is bridged not to
a resonator, but to a real ShortToGround stub (a small vacuum gap away, same
COUPLING_GAP_UM as before) via `create_port_2_conds` with impedance_R=0,
impedance_C=Cc -- i.e. exactly the same lumped-element mechanism, just
grounded instead of routed to a resonator, removing every other confound.
"""
import numpy as np

CPW_WIDTH_UM = 20.0
CPW_GAP_UM = 11.5
Z0_OHM = 49.961
T_LENGTH_UM = 15.0
COUPLING_GAP_UM = 20.0
CHIP_HALF_SPAN_UM = 800.0


def S21_analytic(freq_hz, C_F, Z0=Z0_OHM):
    """Exact S21 for an ideal shunt capacitor C on a matched Z0 line."""
    omega = 2 * np.pi * np.asarray(freq_hz)
    Y = 1j * omega * C_F
    return 1.0 / (1.0 + Y * Z0 / 2.0)


def C_from_S21(freq_hz, S21):
    """Invert the same formula to recover the effective realized C from a
    measured (simulated) S21(f) -- solves Y = 2*(1-S21)/(S21*Z0), C=Y/(j*omega)."""
    omega = 2 * np.pi * np.asarray(freq_hz)
    Y = 2.0 * (1.0 - S21) / (S21 * Z0_OHM)
    return (Y / (1j * omega)).real   # imag part should be ~0 for a pure C; real
                                     # part carries any resistive contamination


def build_design(Cc_F):
    from qiskit_metal import designs, Dict
    from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
    from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
    from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder
    from qiskit_metal.qlibrary.couplers.line_tee import LineTee

    chip_size_x_mm = 2.0 * (CHIP_HALF_SPAN_UM + 400.0) / 1e3
    gnd_y_um = -(T_LENGTH_UM + COUPLING_GAP_UM)
    chip_size_y_mm = (abs(gnd_y_um) + 400.0) / 1e3
    chip_center_y_um = gnd_y_um / 2.0

    design = designs.DesignPlanar({}, overwrite_enabled=True)
    design.chips.main.size.size_x = f'{chip_size_x_mm}mm'
    design.chips.main.size.size_y = f'{chip_size_y_mm}mm'
    design.chips.main.size.size_z = '500um'
    design.chips.main.size.center_x = '0mm'
    design.chips.main.size.center_y = f'{chip_center_y_um*1e-3}mm'
    design.variables['cpw_width'] = f'{CPW_WIDTH_UM}um'
    design.variables['cpw_gap'] = f'{CPW_GAP_UM}um'

    x1 = f'-{CHIP_HALF_SPAN_UM}um'
    x2 = f'{CHIP_HALF_SPAN_UM}um'
    OpenToGround(design, 'OTG1', options=dict(
        chip='main', pos_x=x1, pos_y='0um', orientation='180', width='cpw_width'))
    OpenToGround(design, 'OTG2', options=dict(
        chip='main', pos_x=x2, pos_y='0um', orientation='0', width='cpw_width'))

    LineTee(design, 'JT', options=dict(
        chip='main', pos_x='0um', pos_y='0um', orientation='0',
        prime_width='cpw_width', prime_gap='cpw_gap',
        second_width='cpw_width', second_gap='cpw_gap',
        t_length=f'{T_LENGTH_UM}um'))

    RoutePathfinder(design, 'TL1', options=dict(
        chip='main', trace_width='cpw_width', trace_gap='cpw_gap', fillet='50um',
        pin_inputs=Dict(start_pin=Dict(component='OTG1', pin='open'),
                        end_pin=Dict(component='JT', pin='prime_start'))))
    RoutePathfinder(design, 'TL2', options=dict(
        chip='main', trace_width='cpw_width', trace_gap='cpw_gap', fillet='50um',
        pin_inputs=Dict(start_pin=Dict(component='JT', pin='prime_end'),
                        end_pin=Dict(component='OTG2', pin='open'))))

    ShortToGround(design, 'GND1', options=dict(
        chip='main', pos_x='0um', pos_y=f'{gnd_y_um}um', orientation='90',
        width='cpw_width'))

    design.rebuild()
    return design
