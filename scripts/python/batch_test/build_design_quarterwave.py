"""
Builds a BARE quarter-wave (lambda/4) resonator: open-short topology
(otg1 = OpenToGround, sg1 = ShortToGround), no transmission line or
launchpads - just the resonator, for eigenmode-only studies of a single
grounding strap swept around a target harmonic.

--- Length calibration (5 GHz fundamental) ---
Uses design_helper.DesignHelper's elliptic-integral CPW model directly
for THIS design's actual cross-section (cpw_width=20um, cpw_gap=12.25um),
NOT a cross-topology rescale from a different geometry - that approach
was tried earlier in this project and found to accumulate real error
(a topology rescale from a 10um/6um half-wave measurement was off by
>50% before elliptic-integral correction; see LESSONS_LEARNED).

    eps_r = 11.9          (silicon, raw substrate permittivity - NOT
                            pre-averaged; compute_cpw_geometric_parameters
                            already does the vacuum/substrate averaging
                            internally via the elliptic integrals)
    alpha_inductance = 0  (matches Palace's ideal Dirichlet PEC boundary,
                            no kinetic inductance)
    t = 100nm              (metal thickness; v is only weakly sensitive
                            to this - t=100nm vs 200nm shifts L by <0.5%)

    v = compute_cpw_geometric_parameters('20um', '12.25um', '100nm', 11.9)['v']
      = 1.1858e8 m/s

Open-short quarter-wave: f1 = v / (4L)  =>  L = v / (4 * f1)
    L = 1.1858e8 / (4 * 5e9) = 5.9291 mm

This independently reproduces the 5.92905mm already in this file to
better than 0.03% - the number was right, only the topology (both ends
were OpenToGround, i.e. half-wave, not open-short) needed fixing to
actually match it.

TREAT THIS AS A STARTING ESTIMATE, same as the half-wave design's own
calibration history: fillets, lead-in straights, and the ShortToGround/
OpenToGround termination geometries all introduce second-order effects
the ideal transmission-line formula doesn't capture. Run the ungrounded
baseline (ground_pos_mm=None) FIRST and compare measured f1 against
5.0 GHz before trusting any strap-position sweep on this length - the
half-wave design needed two correction iterations (+54.6% then -3.1%
then -1.7%) before converging; expect something similar here.
"""
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["PMIX_MCA_gds"] = "hash"

from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.meandered_grounded import RouteMeanderGrounded
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround

# Derived length for a 5 GHz quarter-wave (open-short) fundamental at
# THIS design's cpw_width/cpw_gap - see docstring above. Independently
# cross-checked via design_helper.DesignHelper's elliptic-integral CPW
# model to 1.1858e8/(4*5e9)*1000 = 5.9291mm, matching to <0.03%.
TOTAL_LENGTH_MM = 5.92905
GROUND_WIDTH_UM = 10
# How far each strap end extends PAST the outer gap edge into the ground
# plane. With 0 the strap terminates exactly ON the ground-plane boundary,
# a coincident edge that meshes badly - worst on filleted corners, where a
# straight strap rectangle meets a curved gap boundary at a shallow angle
# and leaves a thin wedge. Increase this if sliver elements persist at the
# strap corners; it is cheap to sweep with mesh_test_qw_local.py's
# --strap-overlap-um flag.
GROUND_OVERLAP_UM = 4
# Keep the strap this far from either physical end of the meander path,
# since a strap placed exactly at an endpoint is not physically meaningful.
EDGE_MARGIN_MM = 0.5

# Chip footprint. Deliberately small since this design has no TL/
# launchpads to accommodate - just enough substrate/airbox margin around
# the resonator itself. Re-run the baseline after changing these, same
# caveat as the half-wave design: ground-plane proximity and airbox
# height both measurably shift f1 (the schizo-resonator chip-shrink test
# found a 25% Y-volume cut moved f1 by only 0.045%, but that was a
# DIFFERENT geometry - unverified here until tested).
CHIP_SIZE_X_MM = 3.0
CHIP_SIZE_Y_MM = 1.5
CHIP_SIZE_Z_UM = 500
CHIP_CENTER_X_MM = 0.0
CHIP_CENTER_Y_MM = -0.7


def build_design(ground_pos_mm=None, total_length_mm=TOTAL_LENGTH_MM,
                  chip_size_x_mm=None, chip_size_y_mm=None,
                  chip_center_x_mm=None, chip_center_y_mm=None,
                  ground_overlap_um=None):
    """
    Build the bare quarter-wave resonator design (no TL, no launchpads).

    ground_pos_mm : float or None
        Position (mm) along the resonator path, measured from the SHORT
        end (sg1's pin), for a single grounding strap. If None, builds
        the UNGROUNDED baseline (plain RouteMeander, no strap) - this is
        the reference point to validate/calibrate TOTAL_LENGTH_MM against
        the true measured f1, same role as in the half-wave design.
    total_length_mm : float
        Physical path length of the meander, in mm.
    chip_size_x_mm, chip_size_y_mm, chip_center_x_mm, chip_center_y_mm :
        float or None
        Chip footprint overrides. None uses the CHIP_* module defaults.
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

    # OPEN end (x=L, per the docstring's boundary-condition convention:
    # short at x=0, open at x=L) - this is what makes it lambda/4, NOT
    # OpenToGround at both ends (that would be the half-wave design).
    # width/gap MUST be given explicitly. OpenToGround and ShortToGround
    # carry literal default_options of width='10um', gap='6um' - they do
    # NOT inherit design.variables['cpw_width']/['cpw_gap']. Leaving them
    # unset builds both terminations at the half-wave design's 10/6um
    # while the trace itself is 20/12.25um, producing a visible step at
    # the open end (and a small, silent impedance discontinuity there).
    # termination_gap is the extra gap past the open end that sets the
    # end capacitance - matching it to the CPW gap keeps the open
    # electrically consistent with the line.
    otg1 = OpenToGround(design, 'otg1', options=dict(
        chip='main', pos_x='-0.2mm', pos_y='-40um', orientation=180,
        width='20um', gap='12.25um', termination_gap='12.25um'))
    # SHORT end (x=0). No termination_gap here - the conductor is tied to
    # ground, so there is no open-end gap to size.
    sg1 = ShortToGround(design, 'sg1', options=dict(
        chip='main', pos_x='0mm', pos_y='-1.35mm', orientation=-90,
        width='20um', gap='12.25um'))

    common_kwargs = dict(
        trace_width='20um',
        trace_gap='12.25um',
        total_length=f'{total_length_mm}mm',
        hfss_wire_bonds=False,
        fillet='99.9 um',
        lead=dict(start_straight='300um'),
        pin_inputs=Dict(
            start_pin=Dict(component='sg1', pin='short'),
            end_pin=Dict(component='otg1', pin='open')),
    )

    if ground_pos_mm is None:
        # Baseline: no grounding strap at all - plain quarter-wave resonator
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

    d_grounded = build_design(ground_pos_mm=2.0)
    print("Grounded design built OK. Components:", list(d_grounded.components.keys()))
