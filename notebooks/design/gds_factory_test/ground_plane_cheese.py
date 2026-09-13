# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2017, 2021.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.
"""Triangular-lattice circular-hole ground-plane cheesing, native to QDesign.

Installed into the qcg-quantum-design env's
qiskit_metal/qlibrary/sample_shapes/, alongside this repo's other custom
qlibrary components. This copy here (notebooks/design/gds_factory_test/) is
the demo/dev copy kept in sync with it -- edit the installed copy as the
source of truth if the two drift.
"""

import numpy as np
import shapely
from shapely.ops import unary_union
import klayout.db as kdb

from qiskit_metal import Dict
from qiskit_metal.qlibrary.core import QComponent

_MM_TO_DBU = 1e6  # klayout works in integer nm; qiskit_metal geometry is in mm


def _hex_lattice_centers(width, height, pitch):
    """Centers of a triangular (hex-packed) lattice covering a box of size
    ``width`` x ``height`` centred on the origin, nearest-neighbour spacing
    ``pitch``. All lengths in mm."""
    dx = pitch
    dy = pitch * np.sqrt(3) / 2
    nx = int(np.ceil(width / dx)) + 2
    ny = int(np.ceil(height / dy)) + 2
    centers = []
    for j in range(-ny, ny + 1):
        y = j * dy
        if abs(y) > height / 2 + pitch:
            continue
        xoff = (dx / 2) if (j % 2) else 0.0
        for i in range(-nx, nx + 1):
            x = i * dx + xoff
            if abs(x) <= width / 2 + pitch and abs(y) <= height / 2 + pitch:
                centers.append((x, y))
    return np.array(centers)


def _place_centers(local_centers, pos_x, pos_y, orientation_deg):
    """Rotate ``local_centers`` (Nx2, mm, lattice centred on the origin) by
    ``orientation_deg`` about the origin and translate to (pos_x, pos_y) --
    i.e. the same rigid transform ``draw.rotate``/``draw.translate`` would
    apply to the final hole shapes, but applied to the lattice points
    directly so every downstream step (keepout filtering, hole generation)
    already works in the design's absolute coordinates."""
    theta = np.radians(orientation_deg)
    c, s = np.cos(theta), np.sin(theta)
    rot = np.array([[c, -s], [s, c]])
    return local_centers @ rot.T + np.array([pos_x, pos_y])


def _ngon_hole_region(centers_mm, radius_mm, n_sides, orientation_deg=0.0):
    """klayout Region of one n-gon hole per (absolute-coordinate) center, in
    nm. Holes at a hex-lattice pitch don't overlap, so no merge-induced hole
    shapes appear -- each polygon in the returned region is a single simple
    n-gon. ``orientation_deg`` rotates each hole's own vertices to match the
    lattice's rotation (only visible for a small n_sides)."""
    theta0 = np.radians(orientation_deg)
    angles = 2 * np.pi * np.arange(n_sides) / n_sides + theta0
    template_nm = np.stack(
        [radius_mm * np.cos(angles), radius_mm * np.sin(angles)], axis=1
    ) * _MM_TO_DBU
    region = kdb.Region()
    for cx, cy in centers_mm * _MM_TO_DBU:
        xy = template_nm + np.array([cx, cy])
        region.insert(kdb.Polygon([kdb.Point(int(round(x)), int(round(y))) for x, y in xy]))
    region.merge()
    return region


def _region_to_shapely(region):
    """klayout Region (nm) -> shapely MultiPolygon (mm). Hull points only --
    the hole polygons this module builds are simple n-gons with no interior
    rings, so klayout never hands back a polygon with holes here."""
    polys = [
        shapely.Polygon([(pt.x / _MM_TO_DBU, pt.y / _MM_TO_DBU) for pt in p.each_point_hull()])
        for p in region.each()
    ]
    return shapely.MultiPolygon(polys)


class GroundPlaneCheese(QComponent):
    """Cheese a rectangular patch of the ground plane with a triangular
    (hex-packed) lattice of circular holes, approximated as N-gons to bound
    mesh growth per hole.

    Emits a single ``subtract=True`` poly, same convention as ``RectCutout``:
    every renderer removes it from the ground metal only. Holes whose center
    falls within ``keepout_margin`` (plus their own radius) of ANY existing
    edge already in the design -- metal boundary or another component's own
    cutout boundary (e.g. a ground-plane pocket) -- are dropped whole, never
    clipped into slivers. Because this reads the design's current qgeometry,
    add this component *after* every other component it must keep clear of.

    Cheesing can fragment the ground plane into disconnected islands, which
    breaks the return-current path in a Palace eigenmode solve. This
    component checks that after it's added (via ``qmview.ground_pieces``,
    chip sheet minus every subtract poly in the whole design) and raises
    immediately if the result isn't a single connected piece -- it does not
    silently hand a fragmented plane to the pipeline.

    Note: pairwise-disjoint circular holes (pitch > 2*hole_radius) can never
    disconnect a plane by themselves -- topologically, removing any set of
    disjoint bounded punctures from a connected region leaves it connected.
    This check only fires if holes overlap each other/the outer boundary, or
    combine with other subtract geometry already in the design.

    .. meta::
        :description: Triangular-lattice circular-hole ground-plane cheesing

    Default Options:
        * pos_x, pos_y: '0um' -- center of the cheesed region
        * width, height: '1000um' -- size of the region to cheese
        * pitch: '100um' -- nearest-neighbour hole spacing
        * hole_radius: '20um' -- hole radius (before N-gon approximation)
        * n_sides: '12' -- polygon sides approximating each circle, max 16
        * keepout_margin: '15um' -- clearance kept around existing geometry
        * orientation: '0'
        * helper: 'False'
    """

    default_options = Dict(
        pos_x='0um',
        pos_y='0um',
        width='1000um',
        height='1000um',
        pitch='100um',
        hole_radius='20um',
        n_sides='12',
        keepout_margin='15um',
        orientation='0',
        helper='False',
    )
    """Default options."""

    component_metadata = Dict(short_name='cheese', _qgeometry_table_poly='True')
    """Component metadata."""

    TOOLTIP = "Triangular-lattice circular-hole ground-plane cheesing."

    def make(self):
        """Convert self.options into QGeometry."""
        p = self.parse_options()
        n_sides = int(p.n_sides)
        if not (3 <= n_sides <= 16):
            raise ValueError(
                f"GroundPlaneCheese '{self.name}': n_sides={n_sides} out of "
                "range [3, 16] -- higher side counts defeat the point of "
                "N-gon-approximating the hole (bounding per-hole mesh growth)."
            )

        # Build the lattice centred on the origin, then place it at
        # (pos_x, pos_y)/orientation -- keepout filtering below compares
        # against the design's absolute-coordinate qgeometry, so it MUST run
        # on absolute-coordinate centers, not local ones.
        centers = _hex_lattice_centers(p.width, p.height, p.pitch)
        centers = _place_centers(centers, p.pos_x, p.pos_y, p.orientation)

        # Keep clear of every existing shape -- metal (a real component's
        # island, pad, or CPW trace) and cutout (another component's own
        # ground-plane pocket, or a CPW's gap) need different tests. Traces
        # live in the 'path' qgeometry table as a LineString + width, not as
        # a 'poly' -- read both, or a CPW trace is invisible to this check
        # and a hole can cut straight through its metal with zero clearance.
        #
        # Metal must never be touched anywhere, edge or interior, so it gets
        # a full-polygon buffer: a hole is excluded if it's anywhere within
        # margin+radius of the shape, inside or out. A cutout's *interior* is
        # already void ground -- a hole landing fully inside one is
        # redundant, not wrong -- so it only needs a boundary buffer: excluded
        # only if it would straddle the existing cutout's edge.
        poly = self.design.qgeometry.tables['poly']
        path = self.design.qgeometry.tables['path']
        metal = list(poly[~poly['subtract'].astype(bool)]['geometry'])
        metal += [row.geometry.buffer(row.width / 2)
                  for row in path.itertuples() if not row.subtract]
        cutouts = list(poly[poly['subtract'].astype(bool)]['geometry'])
        cutouts += [row.geometry.buffer(row.width / 2)
                    for row in path.itertuples() if row.subtract]

        keepout_parts = []
        if metal:
            keepout_parts.append(unary_union(metal).buffer(p.keepout_margin + p.hole_radius))
        if cutouts:
            edges = unary_union([g.boundary for g in cutouts])
            keepout_parts.append(edges.buffer(p.keepout_margin + p.hole_radius))
        if keepout_parts:
            keepout = unary_union(keepout_parts)
            centers = centers[~shapely.intersects(shapely.points(centers), keepout)]

        region = _ngon_hole_region(centers, p.hole_radius, n_sides, p.orientation)
        holes = _region_to_shapely(region)

        n_holes = len(holes.geoms)
        region_area = p.width * p.height
        self.cheese_report = dict(
            n_holes=n_holes,
            coverage=holes.area / region_area if region_area else 0.0,
            mesh_element_proxy=n_holes * n_sides,
        )
        print(
            f"GroundPlaneCheese '{self.name}': {n_holes} holes, "
            f"coverage={self.cheese_report['coverage']:.3f}, "
            f"mesh element proxy (n_holes * n_sides)={self.cheese_report['mesh_element_proxy']}"
        )

        self.add_qgeometry(
            'poly', {'cheese_holes': holes}, subtract=True, helper=p.helper,
            layer=p.layer, chip=p.chip,
        )

        import qmview
        pieces = qmview.ground_pieces(self.design)
        if len(pieces) != 1:
            raise RuntimeError(
                f"GroundPlaneCheese '{self.name}': cheesing fragmented the "
                f"ground plane into {len(pieces)} disconnected pieces "
                f"(pitch={p.pitch}mm, hole_radius={p.hole_radius}mm, "
                f"keepout_margin={p.keepout_margin}mm). Widen pitch, shrink "
                "hole_radius, or narrow the cheesed region so bridges survive."
            )
