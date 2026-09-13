"""Convert PalaceForCQED's raw .geo elements 0 (a bare hanger resonator) and
4/5 (the feedline) directly into a GDS file, for feeding into SQDMetal via
gds_design=. Literal 1:1 point conversion -- no re-derivation of dimensions."""
import os
import re
import gdstk

GEO_PATH = "/home/ruiz/Documents/thesis/notebooks/design/parametric_coupling_chip/s11_sqdmetal_validation/palace_cqed_reference/2DQv8_eb4_edit.geo"
OUT_GDS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "phase2_geo_resonator.gds")

with open(GEO_PATH) as f:
    lines = f.readlines()

elem_starts = []
for i, l in enumerate(lines):
    m = re.search(r"element #(\d+)", l)
    if m:
        elem_starts.append((i, int(m.group(1))))
elem_starts.append((len(lines), None))

point_re = re.compile(r"Point\((\d+)\)=\{([-+0-9.eE]+),\s*([-+0-9.eE]+),")

elems = {}
for idx in range(len(elem_starts) - 1):
    start, eid = elem_starts[idx]
    end = elem_starts[idx + 1][0]
    pts = []
    for l in lines[start:end]:
        m = point_re.search(l)
        if m:
            pts.append((float(m.group(2)), float(m.group(3))))
    if pts:
        elems[eid] = pts

# .geo coordinates are in um; GDS library unit defaults to 1um/1nm precision.
lib = gdstk.Library(unit=1e-6, precision=1e-9)
cell = gdstk.Cell("PHASE2_RESONATOR")

METAL_LAYER = 1
for eid in (0, 4, 5):
    pts = elems[eid]
    poly = gdstk.Polygon(pts, layer=METAL_LAYER, datatype=0)
    cell.add(poly)
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    print(f"element {eid}: {len(pts)} pts, bbox x[{min(xs):.1f},{max(xs):.1f}] "
          f"y[{min(ys):.1f},{max(ys):.1f}]")

lib.add(cell)
lib.write_gds(OUT_GDS)
print(f"wrote {OUT_GDS}")
