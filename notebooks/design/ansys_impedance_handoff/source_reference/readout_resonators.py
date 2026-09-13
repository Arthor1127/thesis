"""Importable build_design() for the hanger readout resonator developed in
readout_resonators.ipynb: a straight feedline (LP1 -- TL1 -- CoupledLineTee --
TL2 -- LP2) with a lambda/4 resonator hanging off the CoupledLineTee's
'second_end' pin, routed as a boustrophedon ("snake") down to a ShortToGround
-- see the notebook for the full derivation, the qiskit_metal RouteMeander
pitfalls that motivated the snake (git history), and the resonator-overlap
safety check this reuses verbatim.

Used by run_qc_sim.py (driven S21 sim on the cluster) to build the
(coupling_space_um, coupling_length_um) sweep points -- see gen_qc_sweep_params.py
for the grid and submit_qc_sweep.sh for the array job.
"""
from collections import OrderedDict

import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.tlines.anchored_path import RouteAnchors
from qiskit_metal.qlibrary.tlines.pathfinder import RoutePathfinder
from qiskit_metal.qlibrary.terminations.launchpad_wb import LaunchpadWirebond
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.couplers.coupled_line_tee import CoupledLineTee

# ---- CPW cross-section: chosen for Z0~50 Ohm at cpw_width=20um -- see
# readout_resonators.ipynb's "Verificación de impedancia de 50 Ω" cell
# (DesignHelper.compute_cpw_geometric_parameters: Z0=49.961 Ohm, -0.08% off). ----
CPW_WIDTH_UM = 20.0
CPW_GAP_UM = 11.5

# ---- target resonator length for an 8 GHz lambda/4 fundamental at the above
# cross-section (eps_r=11.45, thickness=100nm) -- computed once in the notebook
# via DesignHelper.compute_cpw_geometric_parameters's phase velocity (NOT
# DesignHelper.compute_resonator_length, which ignores the CPW's eps_eff --
# see the notebook's note on why). Hardcoded here rather than re-imported from
# `parametric_coupling` so this script has no dependency beyond qiskit_metal/
# SQDMetal on the cluster. Recompute if CPW_WIDTH_UM/CPW_GAP_UM change:
#     from parametric_coupling.design_helper import DesignHelper
#     v = DesignHelper().compute_cpw_geometric_parameters(
#         width=CPW_WIDTH_UM, gap=CPW_GAP_UM, thickness=0.1, eps_r=11.45
#     )['v'].to('m/s').magnitude
#     TARGET_RES_LENGTH_MM = v / (4 * 8.0e9) * 1e3
TARGET_RES_LENGTH_MM = 3.77259

FILLET_RADIUS_UM = 50.0          # feedline/coupler fillet
SNAKE_FILLET_UM = 100.0          # resonator snake's own fillet -- raised 40->70->100um
                                 # (still-sharp-looking turns at 70um); stays under half
                                 # of the also-raised SNAKE_ROW_STEP_UM=300um so
                                 # consecutive turn arcs can't overlap (300 - 2*100 =
                                 # 100um of straight left).
CHIP_HALF_SPAN_UM = 1000.0       # LP1/LP2 distance from center (2mm launcher-to-
                                 # launcher). Briefly raised to 6000um (12mm) this
                                 # session to try to detune the bare feedline cavity
                                 # mode away from 8GHz -- reverted after an overnight
                                 # feedline-length x box-size eigenmode sweep
                                 # (submit_feedline_sweep.sh) showed two things: (1)
                                 # the resonator's own mode is extremely stable
                                 # (8.25-8.28GHz, only 30MHz spread) across feedline
                                 # lengths 4-12mm, so length doesn't matter for the
                                 # resonator itself; (2) the "competing" modes near it
                                 # do NOT separate into clean feedline-vs-box families
                                 # (a same-feedline/different-box-size probe showed
                                 # them shifting with BOTH), meaning this no-port
                                 # eigenmode setup's closed PEC boundary is dominating
                                 # the picture with simulation-box artifacts that
                                 # don't exist in the real, port-terminated device.
                                 # Lengthening the feedline was chasing a number this
                                 # tool can't reliably give -- back to the short
                                 # feedline whose real (ported) driven-S21 behavior
                                 # was already characterized (a mild double-dip).
COUPLING_BUDGET_UM = 700.0       # coupling_length_um + down_length_um, held fixed --
                                 # raised from 550um. Traced directly via the qgeometry
                                 # table: CoupledLineTee's own down_length turn (fillet=
                                 # FILLET_RADIUS_UM=50um) renders cleanly rounded at
                                 # down_length_um=250 but SHARP at 150 and 50 -- the real
                                 # safe threshold is well above the naive "2x radius"
                                 # guess. The production sweep's worst case (coupling_
                                 # length_um=400, see gen_qc_sweep_params.py) used to give
                                 # down_length_um=150 (in the broken zone); the new budget
                                 # gives it 300um instead, comfortably above 250.
SNAKE_X_MAX_UM = 1700.0          # raised way up from 500um: with CHIP_HALF_SPAN_UM now
                                 # 6000um, launchpads are far enough away that the whole
                                 # remaining resonator length fits in a SINGLE east-jog-
                                 # west fold (one meander) instead of multiple stacked
                                 # rows -- no need to keep folding sideways just to avoid
                                 # the (now much more distant) launchpads.
SNAKE_ROW_STEP_UM = 300.0        # raised from 200um -- gives the (now single) row-to-row
                                 # jog more vertical room for the larger SNAKE_FILLET_UM.
SNAKE_INITIAL_LEAD_UM = 250.0    # raised from 100um -- same reasoning as SNAKE_MIN_SEG_UM
                                 # below: 100um was shorter than SNAKE_FILLET_UM itself,
                                 # leaving readout_1's own first corner (initial lead into
                                 # the first horizontal leg) unfilleted too.
SNAKE_MIN_SEG_UM = 2.5 * SNAKE_FILLET_UM   # 250um -- the boustrophedon's leftover/tail
                                           # legs (build_snake_anchors' `max(remaining_um,
                                           # min_seg_um)` fallback) used to default to a
                                           # flat 100um, which is now BELOW SNAKE_FILLET_UM
                                           # itself: a leg needs roughly 2x its fillet
                                           # radius of length to have room for a full,
                                           # non-overlapping rounded corner on each end, so
                                           # a too-short leftover leg rendered as a sharp,
                                           # unfilleted turn instead. Tying this to
                                           # SNAKE_FILLET_UM (with a 2.5x margin over the
                                           # strict 2x minimum) keeps that invariant even if
                                           # the fillet radius changes again later, at the
                                           # cost of a small, sweep-wide-constant overshoot
                                           # of the target resonator length whenever a
                                           # leftover leg gets rounded up to this floor.


def _design_bounds_mm(design):
    b = np.array([t.total_bounds for _, t in design.qgeometry.tables.items() if len(t)])
    return (b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max())


def resonator_overlap_um2(design, other_names=("TL1", "TL2", "CLT", "LP1", "LP2")):
    """resonator<->(other_names) overlap area (um^2) -- should always be 0.
    See readout_resonators.ipynb's identical function for the two bugs this
    caught during development (a negative down_length_um, and a jog too small
    for the meander's own sideways reach)."""
    import shapely
    from shapely.ops import unary_union

    def _trace_union(names):
        ids = [design.components[n].id for n in names if n in design.components]
        polys = []
        pt, pp = design.qgeometry.tables["path"], design.qgeometry.tables["poly"]
        m = (~pt["subtract"].astype(bool)) & pt["component"].isin(ids)
        for _, r in pt[m].iterrows():
            w = float(r.get("width") or 0.0) or 0.01
            polys.append(r["geometry"].buffer(w / 2.0, cap_style=2, join_style=2))
        m = (~pp["subtract"].astype(bool)) & pp["component"].isin(ids)
        polys += list(pp[m]["geometry"])
        return unary_union(polys) if polys else shapely.geometry.GeometryCollection()

    res = _trace_union(["readout_1"])
    other = _trace_union(list(other_names))
    if res.is_empty or other.is_empty:
        return 0.0
    return res.intersection(other).area * 1e6


def build_snake_anchors(start_xy_mm, remaining_mm, x_max_um=SNAKE_X_MAX_UM,
                        row_step_um=SNAKE_ROW_STEP_UM, initial_lead_um=SNAKE_INITIAL_LEAD_UM,
                        min_seg_um=SNAKE_MIN_SEG_UM):
    """Boustrophedon waypoints (um, absolute) -- see readout_resonators.ipynb's
    build_snake_anchors for the full derivation. Returns (waypoints_um,
    final_heading), final_heading in {'EAST', 'WEST', 'DOWN'}."""
    x0, y0 = start_xy_mm[0] * 1e3, start_xy_mm[1] * 1e3
    remaining_um = remaining_mm * 1e3 - initial_lead_um
    assert remaining_um > 0, "remaining_mm too short for even the initial lead"

    x, y = x0, y0 - initial_lead_um
    pts = [(x, y)]

    seg = x_max_um if remaining_um >= x_max_um else max(remaining_um, min_seg_um)
    x += seg; remaining_um -= seg
    pts.append((x, y))
    row_dir, heading = +1, 'EAST'

    while remaining_um > 1e-6:
        seg = row_step_um if remaining_um >= row_step_um else max(remaining_um, min_seg_um)
        y -= seg; remaining_um -= seg
        pts.append((x, y))
        heading = 'DOWN'
        if remaining_um <= 1e-6:
            break
        row_dir = -row_dir
        full_seg = 2 * x_max_um
        seg = full_seg if remaining_um >= full_seg else max(remaining_um, min_seg_um)
        x += row_dir * seg; remaining_um -= seg
        pts.append((x, y))
        heading = 'EAST' if row_dir > 0 else 'WEST'

    return pts, heading


def _build_design_once(coupling_space_um, coupling_length_um,
                       target_length_mm=TARGET_RES_LENGTH_MM,
                       coupling_budget_um=COUPLING_BUDGET_UM,
                       chip_size_x_mm=None, chip_size_y_mm=4.0,
                       chip_half_span_um=CHIP_HALF_SPAN_UM):
    """Feedline + hanger resonator at one (coupling_space_um, coupling_length_um)
    point. Raises if the built geometry overlaps itself (see
    resonator_overlap_um2) -- never returns a disfigured design silently.

    chip_half_span_um is the LP1/LP2 distance from centre, i.e. half the
    launcher-to-launcher feedline length. It is a parameter (not just the
    CHIP_HALF_SPAN_UM module constant) so the feedline length can be swept:
    the bare feedline's own open-open cavity harmonics scale as v/(2L), so
    they move with it while the resonator's own mode stays put -- that
    contrast is what identifies which eigenmode is which (see
    submit_feedline_sweep.sh).

    chip_size_x_mm defaults to None -> auto-sized to fit the launchpads
    (each reaches ~780um beyond its pos_x: lead 30 + taper 350 + pad 400)
    plus ~600um of margin to the chip edge for a clean background mesh taper.
    """
    if chip_size_x_mm is None:
        chip_size_x_mm = 2.0 * (chip_half_span_um + 1400.0) / 1000.0
    assert coupling_length_um < coupling_budget_um, (
        f"coupling_length_um={coupling_length_um} must be < coupling_budget_um="
        f"{coupling_budget_um} (down_length_um = budget - coupling_length_um "
        f"must stay positive -- see readout_resonators.ipynb's docstring on why)")
    down_length_um = coupling_budget_um - coupling_length_um

    design = designs.DesignPlanar({}, overwrite_enabled=True)
    design.chips.main.size.size_x = f'{chip_size_x_mm}mm'
    design.chips.main.size.size_y = f'{chip_size_y_mm}mm'
    design.chips.main.size.size_z = '500um'
    design.chips.main.size.center_x = '0mm'
    design.chips.main.size.center_y = '0.0mm'
    design.variables['cpw_width'] = f'{CPW_WIDTH_UM}um'
    design.variables['cpw_gap'] = f'{CPW_GAP_UM}um'

    x1, y1 = f'-{chip_half_span_um}um', '0um'
    x2 = f'{chip_half_span_um}um'
    LaunchpadWirebond(design, 'LP1', options=dict(
        chip='main', pos_x=x1, pos_y=y1, orientation='360',
        lead_length='30um', pad_height='400um', pad_width='250um',
        pad_gap='150um', taper_height='350um'))
    LaunchpadWirebond(design, 'LP2', options=dict(
        chip='main', pos_x=x2, pos_y=y1, orientation='180',
        lead_length='30um', pad_height='400um', pad_width='250um',
        pad_gap='150um', taper_height='350um'))

    CoupledLineTee(design, 'CLT', options=dict(
        pos_x='0um', pos_y='0um', orientation='0', chip='main',
        prime_width='cpw_width', prime_gap='cpw_gap',
        second_width='cpw_width', second_gap='cpw_gap',
        coupling_space=f'{coupling_space_um}um', coupling_length=f'{coupling_length_um}um',
        down_length=f'{down_length_um}um', fillet=f'{FILLET_RADIUS_UM}um', mirror=False,
        open_termination=True))

    RoutePathfinder(design, 'TL1', options=dict(
        chip='main', trace_width='cpw_width', trace_gap='cpw_gap', fillet=f'{FILLET_RADIUS_UM}um',
        hfss_wire_bonds=True, lead=dict(end_straight='0.1mm'),
        pin_inputs=Dict(start_pin=Dict(component='LP1', pin='tie'),
                        end_pin=Dict(component='CLT', pin='prime_start'))))
    RoutePathfinder(design, 'TL2', options=dict(
        chip='main', trace_width='cpw_width', trace_gap='cpw_gap', fillet=f'{FILLET_RADIUS_UM}um',
        hfss_wire_bonds=True, lead=dict(start_straight='0.1mm'),
        pin_inputs=Dict(start_pin=Dict(component='CLT', pin='prime_end'),
                        end_pin=Dict(component='LP2', pin='tie'))))

    design.rebuild()

    consumed_mm = (coupling_length_um + down_length_um) * 1e-3
    remaining_mm = target_length_mm - consumed_mm
    if remaining_mm <= 0.2:
        raise ValueError(f"remaining_mm={remaining_mm:.3f} too short to route "
                         f"(coupling_length_um + down_length_um already consumed "
                         f"{consumed_mm*1e3:.0f}um of the {target_length_mm:.3f}mm target)")

    se = design.components['CLT'].pins['second_end']['middle']
    # Keep the snake clear of the launchpads: its east reach is
    # x_max + coupling_length/2, so cap x_max at 400um inside chip_half_span_um.
    # For the (long) default feedline this is inert -- SNAKE_X_MAX_UM wins and
    # every point folds identically -- it only bites for short-feedline sweep
    # points, which fold into more rows to fit the same total length.
    snake_x_max_um = min(SNAKE_X_MAX_UM, chip_half_span_um - 400.0)
    pts_um, heading = build_snake_anchors(se, remaining_mm, x_max_um=snake_x_max_um)

    orientation = {'EAST': '0', 'WEST': '180', 'DOWN': '270'}[heading]
    end_x_um, end_y_um = pts_um[-1]
    ShortToGround(design, 'res_short', options=dict(
        pos_x=f'{end_x_um}um', pos_y=f'{end_y_um}um', orientation=orientation,
        chip='main', width='cpw_width'))

    anchors = OrderedDict()
    for i, (x_um, y_um) in enumerate(pts_um[:-1]):
        anchors[i] = np.array([x_um, y_um]) * 1e-3

    RouteAnchors(design, 'readout_1', options=dict(
        pin_inputs=dict(
            start_pin=dict(component='CLT', pin='second_end'),
            end_pin=dict(component='res_short', pin='short'),
        ),
        fillet=f'{SNAKE_FILLET_UM}um',
        anchors=anchors,
        lead=dict(start_straight='2um', end_straight='2um'),
        trace_width='cpw_width', trace_gap='cpw_gap',
    ))

    design.rebuild()
    overlap = resonator_overlap_um2(design)
    if overlap > 1e-6:
        raise RuntimeError(f"build_design: readout_1 overlaps the coupler/TL/launchpads "
                           f"by {overlap:.2f}um^2 at coupling_space_um={coupling_space_um} "
                           f"coupling_length_um={coupling_length_um}.")

    return design


def build_design(coupling_space_um, coupling_length_um, *args, **kwargs):
    """Thin retry wrapper around _build_design_once. qiskit_metal's RouteAnchors
    has a known fragility (its internal fillet/obstruction-check geometry,
    anchored_path.py's unobstructed_close_up, dies on a None boundary geometry
    for an exact floating-point coincidence in CoupledLineTee's own polygon at
    a sparse, isolated set of (coupling_space_um, coupling_length_um) pairs --
    empirically 1/35 of this project's actual sweep grid, e.g. (10, 350)).
    Nudging readout_1's own anchors by up to 20um did NOT dodge it, pointing at
    the degeneracy being in CLT's geometry rather than readout_1's -- so this
    retries the WHOLE build with coupling_length_um perturbed by a few hundred
    nm instead, which changes CLT's own polygon and reliably avoids the exact
    coincidence, with no physically meaningful effect on Qc/f_r.
    """
    for nudge_um in (0.0, 0.01, 0.1, 0.5, 1.0):
        try:
            return _build_design_once(coupling_space_um, coupling_length_um + nudge_um,
                                      *args, **kwargs)
        except AttributeError as e:
            if "'NoneType' object has no attribute 'coords'" not in str(e):
                raise
    raise RuntimeError(f"build_design: RouteAnchors kept hitting the known "
                       f"None-boundary fragility at coupling_space_um={coupling_space_um} "
                       f"coupling_length_um={coupling_length_um} even after retrying with "
                       f"coupling_length_um nudged by up to 1nm.")


if __name__ == "__main__":
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    for cs, cl in [(10.0, 150.0), (15.0, 200.0), (30.0, 400.0)]:
        d = build_design(cs, cl)
        print(f"cs={cs}um cl={cl}um -> OK, components={list(d.components.keys())}, "
              f"chip={d.chips.main.size.size_x} x {d.chips.main.size.size_y}")
