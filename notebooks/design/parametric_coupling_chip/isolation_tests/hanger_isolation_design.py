"""Geometry for the two-port hanger-topology Qc isolation test.

WHY THIS EXISTS: neither prior isolation test exercises what the real
production sweep actually needs. `run_bareres_eigenmode_isolation.py` is a
1-port DIRECTLY-loaded resonator (analytic Qc via a series R) -- it validated
the eigenmode extraction's numerical precision as a function of Qc, but has
only one lossy port, so it says nothing about SQDMetal/Palace's TWO-PORT
`port-Q.csv` kappa-summing path, which is exactly what showed a ~60x
loaded_Q-vs-kappa disagreement on the real (messy) design. This design fixes
that: a real 2-port hanger (through-line + side-coupled resonator), but with
the coupling realized as an IDEALIZED LUMPED CAPACITOR of a known, exact
value instead of the real design's geometric CoupledLineTee (whose effective
Cc is exactly what's uncertain there) -- so there is still a clean, closed-
form analytic target, this time for the actual 2-port topology.

METHOD: a `LineTee` (a plain galvanic "+"-shaped junction -- NOT
CoupledLineTee, no near-field coupling model at all) gives a real, physically
continuous tap point on the through-line via its `second_end` pin, with a
short (t_length, kept small so its own series inductance is negligible at
8GHz) galvanic stub. The resonator sits a small, real vacuum gap away from
that stub's end, and `create_port_2_conds` bridges the two with an IDEALIZED
lumped capacitor (impedance_C=Cc, impedance_R=0) -- this is the SAME
mechanism used elsewhere in this repo for Josephson-junction lumped
inductors, just with C instead of L, and it draws its own synthetic port
rectangle between the two pins regardless of the real gap geometry, so the
gap size itself has no effect on the applied Cc (same principle as a lumped
port's R being independent of its physical footprint).

ANALYTIC TARGET: standard hanger/notch-filter formula (Pozar; Goppl et al.,
J. Appl. Phys. 104, 113904 (2008); same formula family as this repo's own
design_helper.compute_resonator_Qc_from_Cc, generalized here from a
capacitive-only load's perspective):
    C_eff = pi / (4 * omega_r * Z0)                  (quarter-wave resonator)
    G_ext = (omega_r * Cc)^2 * (Z0 / 2)               (coupling to a MATCHED
                                                        through-line -- Z0/2
                                                        because the matched
                                                        line looks like Z0 in
                                                        EACH direction from the
                                                        tap point, in parallel)
    Qc = omega_r * C_eff / G_ext = 2*C_eff / (Z0 * omega_r * Cc^2)

Three Cc values are pre-solved for Qc_analytic = 300, 3000, 30000 -- a ladder
chosen so the WEAKEST coupling (highest Qc) still gives a comfortably
resolvable driven-S21 notch (linewidth fr/Qc = 267kHz at Qc=30000, still wide
enough for PROM/mesh to resolve, unlike the real design's ~1kHz-scale
effective linewidths at Qc~1e4-1e5), while still climbing into the same
regime where the 1-port isolation test found the eigenmode method's
precision starting to degrade.
"""
import numpy as np

CPW_WIDTH_UM = 20.0
CPW_GAP_UM = 11.5
Z0_OHM = 49.961
F0_HZ = 8.0e9
TARGET_RES_LENGTH_MM = 3.77259   # same lambda/4-at-8GHz length as elsewhere in this repo
T_LENGTH_UM = 15.0               # LineTee stub length -- kept small (<< guided
                                 # wavelength ~14mm at 8GHz) so its own series
                                 # inductance/phase is a negligible perturbation
COUPLING_GAP_UM = 20.0           # real vacuum gap bridged by the lumped-C port
CHIP_HALF_SPAN_UM = 800.0        # OTG1/OTG2 distance from the junction

C_EFF_F = np.pi / (4 * 2 * np.pi * F0_HZ * Z0_OHM)


def Cc_for_Qc(Qc_target):
    """Solve the hanger formula for the coupling capacitance (F) needed to
    hit a target Qc, at this module's Z0/F0."""
    w0 = 2 * np.pi * F0_HZ
    return np.sqrt(2 * C_EFF_F / (Z0_OHM * w0 * Qc_target))


def build_design(Cc_F):
    """Two-port hanger isolation design: OTG1--TL1--[LineTee]--TL2--OTG2
    through-line, with a resonator (OTGr--RES--STG) capacitively tapped off
    the LineTee's `second_end` pin via an idealized lumped capacitor Cc_F.
    Returns (design, junction_component_name, resonator_component_name).
    """
    from qiskit_metal import designs, Dict
    from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
    from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
    from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder
    from qiskit_metal.qlibrary.couplers.line_tee import LineTee

    chip_size_x_mm = 2.0 * (CHIP_HALF_SPAN_UM + 400.0) / 1e3
    res_start_y_um = -(T_LENGTH_UM + COUPLING_GAP_UM)
    chip_size_y_mm = (abs(res_start_y_um) + TARGET_RES_LENGTH_MM * 1e3 + 400.0) / 1e3
    chip_center_y_um = res_start_y_um / 2.0 - TARGET_RES_LENGTH_MM * 1e3 / 2.0

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

    OpenToGround(design, 'OTGr', options=dict(
        chip='main', pos_x='0um', pos_y=f'{res_start_y_um}um', orientation='90',
        width='cpw_width'))
    ShortToGround(design, 'STG', options=dict(
        chip='main', pos_x='0um', pos_y=f'{res_start_y_um - TARGET_RES_LENGTH_MM*1e3}um',
        orientation='270', width='cpw_width'))
    RoutePathfinder(design, 'RES', options=dict(
        chip='main', trace_width='cpw_width', trace_gap='cpw_gap', fillet='50um',
        pin_inputs=Dict(start_pin=Dict(component='OTGr', pin='open'),
                        end_pin=Dict(component='STG', pin='short'))))

    design.rebuild()
    return design
