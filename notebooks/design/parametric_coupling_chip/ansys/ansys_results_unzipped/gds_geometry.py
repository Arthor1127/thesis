"""Read clean metal polygons out of the handoff GDS files.

The handoff GDS files carry the real metal on GDS layer 1, datatypes:
    0  -> ground plane (solid, CPW slot only; NOT cheesed)
    10 -> launchpad pad/taper metal
    11 -> CPW centre-trace metal
Datatypes 99/100/101/102 are cheesing / no-cheese / hole-prototype artifacts
and are ignored -- reading only {0, 10, 11} yields clean solid metal, so the
(heavily cheesed) GDS never has to be imported into HFSS as-is.

qiskit_metal's GDS export stores coordinates whose numeric values are in
millimetres (chip 4.8 -> 4.8 mm, trace half-width 0.01 -> 10 um), despite the
file's unit tag. This module scales by 1000 so everything it returns is in
microns.
"""
import numpy as np
import gdspy

DT_GROUND = 0
DT_PADS = 10
DT_TRACE = 11
_KEEP = (DT_GROUND, DT_PADS, DT_TRACE)

MM_TO_UM = 1000.0  # gdspy returns qiskit_metal coords whose values are in mm


def read_metal_polygons(gds_path, top_cell="TOP_main_1"):
    """Return dict: {'ground': [Nx2 arrays], 'pads': [...], 'trace': [...]} in um."""
    lib = gdspy.GdsLibrary(infile=gds_path)
    if top_cell in lib.cells:
        cell = lib.cells[top_cell]
    else:
        # fall back to the library top level
        tops = lib.top_level()
        cell = tops[0]
    by_spec = cell.get_polygons(by_spec=True)
    out = {"ground": [], "pads": [], "trace": []}
    name = {DT_GROUND: "ground", DT_PADS: "pads", DT_TRACE: "trace"}
    for (layer, dt), polys in by_spec.items():
        if dt not in _KEEP:
            continue
        for p in polys:
            out[name[dt]].append(np.asarray(p, dtype=float) * MM_TO_UM)
    gdspy.current_library = gdspy.GdsLibrary()  # reset global state for repeat calls
    return out


def heal(polys_um):
    """Union a list of Nx2 point arrays into a valid shapely geometry, fixing
    the 'keyhole' polygons qiskit_metal's GDS ground plane uses (a single closed
    ring that dives in to trace each CPW slot and back out along the same line).
    `buffer(0)` turns those into proper polygons-with-holes."""
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    if not polys_um:
        return None
    return unary_union([Polygon(np.asarray(p, float)).buffer(0) for p in polys_um])


def as_polygon_list(geom):
    """Flatten a shapely geometry to a list of Polygons (each may have holes)."""
    if geom is None or geom.is_empty:
        return []
    if geom.geom_type == "Polygon":
        return [geom]
    if geom.geom_type == "MultiPolygon":
        return list(geom.geoms)
    return [g for g in getattr(geom, "geoms", []) if g.geom_type == "Polygon"]


def bbox(polys):
    if not polys:
        return None
    allpts = np.vstack(polys)
    return (allpts[:, 0].min(), allpts[:, 1].min(),
            allpts[:, 0].max(), allpts[:, 1].max())


def summarize(gds_path):
    p = read_metal_polygons(gds_path)
    print(f"--- {gds_path}")
    for k in ("ground", "pads", "trace"):
        b = bbox(p[k])
        n = len(p[k])
        if b is None:
            print(f"  {k:7s}: (none)")
        else:
            print(f"  {k:7s}: {n:3d} poly(s)  x[{b[0]:.1f},{b[2]:.1f}] "
                  f"y[{b[1]:.1f},{b[3]:.1f}] um")
    return p


if __name__ == "__main__":
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    geo = os.path.join(here, "..", "geometry")
    for f in ("line_only.gds", "with_launchpads.gds", "full_design_cs14_cl300.gds"):
        summarize(os.path.join(geo, f))
