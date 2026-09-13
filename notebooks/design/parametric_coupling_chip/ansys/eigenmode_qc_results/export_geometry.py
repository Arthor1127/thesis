"""Export the two new hanger-resonator geometries requested in
../EIGENMODE_QC_REQUEST.md (Part 2) to GDS, using the same generator that
produced ../geometry/full_design_cs14_cl300.gds:
`source_reference/readout_resonators.py::build_design(coupling_space_um,
coupling_length_um)`.

Geometry-only (qiskit_metal, no PALACE / no HFSS) -- safe to run locally.

  conda activate qiskit_qc
  python export_geometry.py                 # writes both new .gds
  python export_geometry.py --check         # also re-export cs14_cl300 and
                                            #   bbox-diff it against the shipped file

Writes: ../geometry/full_design_cs22p5_cl250.gds
        ../geometry/full_design_cs25_cl100.gds
"""
import argparse
import os
import sys

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "source_reference"))
sys.path.insert(0, os.path.join(ROOT, "ansys_results"))

from readout_resonators import build_design          # noqa: E402
from gds_geometry import read_metal_polygons, bbox    # noqa: E402

GEO = os.path.join(ROOT, "geometry")

# (coupling_space_um, coupling_length_um, output basename)
POINTS = [
    (22.5, 250.0, "full_design_cs22p5_cl250"),
    (25.0, 100.0, "full_design_cs25_cl100"),
]


def export(cs_um, cl_um, name):
    design = build_design(cs_um, cl_um)
    out = os.path.join(GEO, name + ".gds")
    design.renderers.gds.export_to_gds(out)
    p = read_metal_polygons(out)
    print(f"wrote {out}")
    for k in ("ground", "pads", "trace"):
        b = bbox(p[k])
        print(f"   {k:7s}: {len(p[k]):2d} poly  "
              + ("x[%.1f,%.1f] y[%.1f,%.1f] um" % (b[0], b[2], b[1], b[3])
                 if b else "(none)"))
    return out


def check_cs14():
    """Re-export cs14_cl300 and bbox-diff against the shipped GDS."""
    shipped = os.path.join(GEO, "full_design_cs14_cl300.gds")
    tmp = os.path.join(HERE, "_cs14_cl300_regen.gds")
    design = build_design(14.0, 300.0)
    design.renderers.gds.export_to_gds(tmp)
    a = read_metal_polygons(shipped)
    b = read_metal_polygons(tmp)
    print("bbox diff (shipped vs regenerated) full_design_cs14_cl300:")
    ok = True
    for k in ("ground", "pads", "trace"):
        ba, bb = bbox(a[k]), bbox(b[k])
        d = max(abs(x - y) for x, y in zip(ba, bb)) if (ba and bb) else None
        print(f"   {k:7s}: shipped {ba}  regen {bb}  max|d|={d}")
        if d is None or d > 1.0:
            ok = False
    os.remove(tmp)
    print("MATCH" if ok else "MISMATCH -- generator differs from the shipped file")
    return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        check_cs14()
    for cs, cl, name in POINTS:
        export(cs, cl, name)
