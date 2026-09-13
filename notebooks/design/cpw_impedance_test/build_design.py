"""A single straight CPW segment (OpenToGround -- RouteStraight -- ShortToGround)
used to validate ``parametric_coupling.design_helper.DesignHelper``'s analytic CPW
formula (impedance Z0, propagation velocity v) against SQDMetal/Palace FEM sims.

Method: run a capacitance sim (electrostatic Maxwell C matrix, trace floating -- both
ends OpenToGround) and a self-inductance sim (magnetostatic flux linkage, current
driven at the open end, returned to ground at a real ShortToGround) on two variants of
the same straight-CPW geometry -- see ``build_design``'s ``far_end`` parameter, which
is NOT optional/cosmetic: get it backwards and the capacitance sim measures the
capacitance of a single shorted conductor (undefined/degenerate), not the CPW line.
Divide each result by the trace length to get C' and L' per unit length, then

    Z0 = sqrt(L'/C')      v = 1/sqrt(L'*C')

exactly the two derived quantities ``compute_cpw_geometric_parameters`` predicts from
(cpw_width, cpw_gap, thickness, eps_r) alone -- see ``analytic_reference`` below.

The metal is the SQDMetal default zero-thickness ("quasi-3D") sheet -- fast and exact
for the in-plane field, matching how every other capacitance sim in this repo treats
thin-film metal (see ``scripts/python/island_cap/run_cap_island.py``). This means the
thickness-dependent finite-t correction terms in Eqs. (B6)-(B7) are not exercised by
the FEM side of this comparison; ``analytic_reference`` still needs *some* thickness
to evaluate the formula (it cannot take t=0), so it defaults to a thin representative
film (100 nm) and the resulting few-% t-dependent term should be treated as part of
the residual between sim and formula, not a bug.

Caveat (fringing / end effects): dividing a FINITE segment's total C and L by its
length is only exact in the limit length >> (cpw_width + 2*cpw_gap); at finite length
each carries an O(1) end correction (open/short termination fringe for C, port/short
fringe for L) that does not shrink with length in absolute terms, only relative to it.
The default LENGTH_UM is set large compared to the swept width/gap range to keep this
a small systematic; ``integration_area_gap_strips`` additionally trims MARGIN_UM off
each end of the flux integral specifically to keep it clear of the port/short field
singularities (a meshing/well-posedness precaution, not primarily an accuracy fix --
the capacitance side has no such lever since Palace only returns one lumped C).
"""
import numpy as np
from qiskit_metal import designs, Dict
from qiskit_metal.qlibrary.terminations.open_to_ground import OpenToGround
from qiskit_metal.qlibrary.terminations.short_to_ground import ShortToGround
from qiskit_metal.qlibrary.tlines.straight_path import RouteStraight

# ---- default geometry ----
LENGTH_UM = 3000.0          # trace length: >> max(cpw_width + 2*cpw_gap) in the sweep
CHIP_MARGIN_UM = 400.0      # airbox margin around the drawn geometry, each side
CHIP_SIZE_Z_UM = 500.0      # substrate thickness

# default trace start/end names, reused by the run_*.py sims
OTG_NAME, STG_NAME, ROUTE_NAME = 'otg', 'stg', 'cpw'

# analytic reference defaults (see module docstring)
REF_THICKNESS_NM = 100.0
REF_EPS_R = 11.45            # SQDMetal Material('silicon'): eps_r = 11.45


def _design_bounds_mm(design):
    b = np.array([t.total_bounds for _, t in design.qgeometry.tables.items() if len(t)])
    return (b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max())


def _size_chip(design, margin_um=CHIP_MARGIN_UM, size_z_um=CHIP_SIZE_Z_UM):
    minx, miny, maxx, maxy = _design_bounds_mm(design)
    m = margin_um * 1e-3
    cx, cy = 0.5 * (minx + maxx), 0.5 * (miny + maxy)
    design.chips.main.size.size_x = f'{(maxx - minx) + 2 * m}mm'
    design.chips.main.size.size_y = f'{(maxy - miny) + 2 * m}mm'
    design.chips.main.size.size_z = f'{size_z_um}um'
    design.chips.main.size.center_x = f'{cx}mm'
    design.chips.main.size.center_y = f'{cy}mm'


def build_design(cpw_width_um, cpw_gap_um, length_um=LENGTH_UM, far_end='short',
                 chip_margin_um=CHIP_MARGIN_UM, chip_size_z_um=CHIP_SIZE_Z_UM):
    """A straight CPW trace of the given width/gap and length, open at x=0.

    ``far_end`` controls the termination at x=length_um, and MUST match the sim
    that will use this design:

    - ``'short'`` (default) -- ``ShortToGround``: a REAL galvanic short (the
      trace polygon runs flush into the ground plane, no isolating gap poly is
      drawn there). Required by run_inductance.py: the self-inductance method
      needs a genuine closed current loop (driven at the open end via a
      U-clip, returned through ground at this real short).
    - ``'open'`` -- ``OpenToGround`` at both ends instead: the trace is a
      single floating conductor, galvanically isolated from ground everywhere.
      Required by run_capacitance.py -- with ``far_end='short'`` the trace and
      ground plane are ONE electrically connected conductor and there is no
      trace<->ground capacitance to extract (Palace would see a single PEC
      region, not two terminals). This was caught after the first version of
      this script shipped with 'short' hardcoded for both sims -- see git
      history / the session that added this docstring paragraph.

    Both terminations use ``termination_gap = cpw_gap`` (>= ~15-20um is the
    safe floor for GMSH to resolve it -- see the
    ``sqdmetal-termination-gap-meshing`` lesson; if you sweep cpw_gap below
    that, fine-mesh the terminations explicitly in the sim script, which
    run_capacitance.py / run_inductance.py already do).
    """
    d = designs.DesignPlanar({}, overwrite_enabled=True)
    d.variables['cpw_width'] = f'{cpw_width_um}um'
    d.variables['cpw_gap'] = f'{cpw_gap_um}um'

    OpenToGround(d, OTG_NAME, options=Dict(
        pos_x='0um', pos_y='0um', orientation='180', chip='main',
        width='cpw_width', gap='cpw_gap', termination_gap='cpw_gap'))

    if far_end == 'short':
        ShortToGround(d, STG_NAME, options=Dict(
            pos_x=f'{length_um}um', pos_y='0um', orientation='0', chip='main',
            width='cpw_width'))
        far_pin = 'short'
    elif far_end == 'open':
        OpenToGround(d, STG_NAME, options=Dict(
            pos_x=f'{length_um}um', pos_y='0um', orientation='0', chip='main',
            width='cpw_width', gap='cpw_gap', termination_gap='cpw_gap'))
        far_pin = 'open'
    else:
        raise ValueError(f"far_end must be 'short' or 'open', got {far_end!r}")

    RouteStraight(d, ROUTE_NAME, options=Dict(
        pin_inputs=Dict(
            start_pin=Dict(component=OTG_NAME, pin='open'),
            end_pin=Dict(component=STG_NAME, pin=far_pin),
        ),
        trace_width='cpw_width', trace_gap='cpw_gap',
        lead=Dict(start_straight='0um', end_straight='0um'),
    ))

    d.rebuild()
    _size_chip(d, chip_margin_um, chip_size_z_um)
    return d


def integration_area_gap_strips(cpw_width_um, cpw_gap_um, length_um=LENGTH_UM,
                                margin_um=None, grow_um=0.3):
    """The two CPW gap strips (trace<->ground, one each side) as a shapely
    MultiPolygon in METRES, for use as a PALACE_Inductance_Simulation integration
    area. ``margin_um`` is trimmed off both x-ends of ``length_um`` to keep the
    integral clear of the open/short termination fields (default: max(5*cpw_gap,
    30um)). ``grow_um`` pushes each strip slightly past the trace/ground edges
    (into the PEC, where B~0) so the integration-area boundary never exactly
    coincides with a metal edge -- an exact coincidence makes GMSH emit degenerate
    zero-length edges when it partitions the surface (see
    ``build_design.loop_aperture_polys`` in flux_optimizing for the same trick).
    """
    import shapely

    if margin_um is None:
        margin_um = max(5.0 * cpw_gap_um, 30.0)
    if 2 * margin_um >= length_um:
        raise ValueError(f"margin_um={margin_um} trims away the whole "
                         f"length_um={length_um} segment; lengthen the trace or "
                         f"shrink margin_um.")

    x0, x1 = margin_um * 1e-6, (length_um - margin_um) * 1e-6
    w, g, gr = cpw_width_um * 1e-6, cpw_gap_um * 1e-6, grow_um * 1e-6
    y_in, y_out = w / 2.0 - gr, w / 2.0 + g + gr

    top = shapely.geometry.box(x0, y_in, x1, y_out)
    bot = shapely.geometry.box(x0, -y_out, x1, -y_in)
    return shapely.MultiPolygon([top, bot])


def analytic_reference(cpw_width_um, cpw_gap_um, thickness_nm=REF_THICKNESS_NM,
                       eps_r=REF_EPS_R):
    """DesignHelper's analytic Cg, Lg, Z0, v, eps_eff for this geometry -- the
    prediction the sim results are validated against. Local import: keeps
    build_design's geometry-building path free of a hard dependency on the
    `parametric-coupling` package, which run_capacitance.py / run_inductance.py
    (executed on the cluster) do not need.
    """
    from parametric_coupling.design_helper import DesignHelper

    dh = DesignHelper()
    return dh.compute_cpw_geometric_parameters(
        width=cpw_width_um, gap=cpw_gap_um, thickness=thickness_nm * 1e-3,
        eps_r=eps_r)


if __name__ == "__main__":
    d = build_design(cpw_width_um=20.0, cpw_gap_um=12.0, length_um=1000.0)
    print(f"components={list(d.components.keys())}  "
          f"chip={d.chips.main.size.size_x} x {d.chips.main.size.size_y}")
    ref = analytic_reference(20.0, 12.0)
    print(f"analytic @ w=20um g=12um: Z0={ref['Z0']:.3f}  v={ref['v'].to('m/s'):.4e}  "
          f"eps_eff={ref['eps_eff']:.3f}")
