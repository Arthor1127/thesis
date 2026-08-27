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

import numpy as np
from qiskit_metal import Dict, draw
from qiskit_metal.qlibrary.tlines.meandered import RouteMeander
from qiskit_metal.renderers.renderer_mpl.mpl_renderer import QMplRenderer


class RouteMeanderGrounded(RouteMeander):
    """Same as `RouteMeander`, but with a set of rectangular grounding straps
    placed along the length of the resonator, each shorting the center
    conductor to the ground plane. Useful e.g. to build periodically-grounded
    resonators.

    Inherits `RouteMeander` class.

    Default Options:
        * ground_straps: Dict
            * positions: '' -- List of lengths (e.g. ``['1mm', '2.5mm']``),
              measured as arc length along the resonator starting from its
              `start_pin`, at which a grounding strap is placed.
            * overlap: '4um' -- How far each end of the strap extends PAST the
              outer edge of the CPW gap, into the ground plane. Must be > 0.
              With overlap=0 the strap terminates exactly ON the ground-plane
              boundary, a coincident-edge condition: on a straight run that is
              a degenerate zero-width sliver, and on a filleted corner it is
              worse, because the strap is a straight rectangle while the gap
              boundary is an arc, so the two meet at a shallow angle and leave
              a thin wedge at each corner. Those wedges force GMSH to generate
              sliver elements. A few um of overlap removes the coincidence
              entirely. It also covers the straight-rectangle-vs-arc error,
              which is only ~w^2/(8*R) (about 0.13um for a 10um strap on a
              100um fillet) - far too small to justify building the strap as
              a curved polygon.
            * width: '10um' -- Width of each strap, measured along the path
              direction (i.e. how "thick" the short is along the resonator).
              Either a single value (used for every strap) or a list matching
              `positions` one-to-one, to give each strap its own width.
    """

    component_metadata = Dict(short_name='cpw_grounded')
    """Component metadata"""

    default_options = Dict(ground_straps=Dict(positions='', width='10um',
                                               overlap='4um'))
    """Default options"""

    TOOLTIP = """Implements a RouteMeander with periodic grounding straps along its length."""

    def make(self):
        """Builds the meander the same way `RouteMeander` does, then adds the
        grounding straps on top."""
        super().make()
        self.make_ground_straps(self.get_points())

    def make_ground_straps(self, pts: np.ndarray):
        """Adds rectangular metal straps shorting the trace to the ground
        plane, at the arc-length positions given in
        `options.ground_straps.positions`.

        Args:
            pts (np.ndarray): Route points, the same array used to build the trace path.
        """
        positions = self.p.ground_straps.positions
        if positions is None or len(positions) == 0:
            return

        widths = self.p.ground_straps.width
        if not isinstance(widths, (list, np.ndarray)):
            widths = [widths] * len(positions)
        if len(widths) != len(positions):
            raise ValueError(
                "options.ground_straps.width must be either a single value "
                "or a list of the same length as options.ground_straps.positions"
            )

        line = self._get_filleted_centerline(pts)
        # Spans the trace and both gaps, fully shorting the centre conductor
        # to both ground planes, and then extends PAST both gap edges into
        # the ground plane rather than
        # stopping exactly on the boundary - see 'overlap' in the class
        # docstring for why the exact-fit case meshes badly, especially
        # on filleted corners.
        overlap = self.p.ground_straps.overlap
        strap_length = self.p.trace_width + 2 * self.p.trace_gap + 2 * overlap

        straps = {}
        for i, (position, width) in enumerate(zip(positions, widths)):
            if not 0 <= position <= line.length:
                raise ValueError(
                    f"ground strap position {position} falls outside of the "
                    f"resonator's length (0, {line.length})")
            straps[f'ground_strap_{i}'] = self._make_strap(
                line, position, width, strap_length)

        self.add_qgeometry('poly', straps)

    def _get_filleted_centerline(self, pts: np.ndarray) -> draw.LineString:
        """Returns the resonator's actual rendered centerline, i.e. `pts` with
        rounded corners applied (radius `self.p.fillet`), same as what
        `make_elements` draws and what SQDMetal meshes.

        With a fillet radius comparable to (or larger than) the meander
        spacing - as is typical for these resonators - most of the visible
        curve sits inside a rounded corner, so interpolating on the raw,
        un-rounded `pts` polyline places straps measurably off the real
        trace and with the wrong tangent direction.

        Args:
            pts (np.ndarray): Route points (straight-segment polyline, pre-fillet).

        Returns:
            shapely.geometry.LineString: The rounded-corner centerline.
        """
        raw_line = draw.LineString(pts)
        if self.p.fillet <= 0:
            return raw_line
        mpl_renderer = QMplRenderer(design=self.design, canvas=None, logger=None)
        return mpl_renderer.fillet_path({'geometry': raw_line, 'fillet': self.p.fillet})

    @staticmethod
    def _make_strap(line: draw.LineString, position: float, width: float,
                    strap_length: float) -> draw.Polygon:
        """Builds a single rectangle crossing the CPW gap at the given
        arc-length position along `line`, oriented perpendicular to the
        path's local direction.

        Args:
            line (LineString): Centerline of the resonator trace.
            position (float): Arc-length distance along `line` at which to place the strap.
            width (float): Strap width, measured along the path direction.
            strap_length (float): Strap length, measured across the trace and both gaps.

        Returns:
            shapely.geometry.Polygon: The strap rectangle, placed and rotated in the design frame.
        """
        # local tangent direction, estimated from two nearby points on the line
        eps = max(line.length * 1e-6, 1e-9)
        p_before = line.interpolate(max(position - eps, 0))
        p_after = line.interpolate(min(position + eps, line.length))
        tangent = np.array([p_after.x - p_before.x, p_after.y - p_before.y])
        tangent = tangent / np.linalg.norm(tangent)
        angle = np.degrees(np.arctan2(tangent[1], tangent[0]))

        center = line.interpolate(position)

        # rectangle(width, strap_length) is built with `width` along local x and
        # `strap_length` along local y; rotating by `angle` (the tangent's angle
        # from +x) aligns local-x with the path direction and local-y with the
        # perpendicular (normal) direction, so the strap crosses the gaps sideways
        rect = draw.rectangle(width, strap_length)
        rect = draw.rotate(rect, angle)
        rect = draw.translate(rect, center.x, center.y)
        return rect

