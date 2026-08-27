"""
Step 1 of the field-visualization workflow (see field_viz_by_component.py
for step 2). Prints a design component's bounding box in MILLIMETERS (matching the
.pvd file's actual coordinate units - see get_component_bbox_mm() below
for why this is mm and not meters), ready to
copy-paste as --xmin/--ymin/--xmax/--ymax into the pvpython script.

Runs under your normal conda env (qcg-quantum-design) since it needs
qiskit_metal/SQDMetal - pvpython does NOT have those installed, which is
why this is a separate script rather than one combined tool.

Usage:
    conda activate qcg-quantum-design
    python get_component_bbox.py --ground-pos 5.906 --component resonator1
    python get_component_bbox.py --ground-pos 5.906 --list-components
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from build_design import build_design, TOTAL_LENGTH_MM


def parse_ground_pos(s):
    if s.lower() == "none":
        return None
    return float(s)


def get_component_bbox_mm(design, component_name):
    """Returns (xmin,ymin,xmax,ymax) in MILLIMETERS for the given
    component. NOTE: despite Palace's own console log reporting mesh
    bounds in meters (e.g. 'Xmin = -2.880e-03 m'), the actual .pvd/.vtu
    file it writes stores raw coordinates in the SAME units GMSH/
    qiskit-metal used to build the mesh - millimeters - NOT meters.
    Confirmed by comparing ParaView's GetBounds() (X=[-2.88,2.88]) against
    the Palace log's Zmin/Zmax*1000 for the same run. Do NOT convert to
    meters here or the clip box will be 1000x too small and silently grab
    the wrong region near the origin instead of erroring out."""
    if component_name == "ground_strap":
        poly_table = design.qgeometry.tables["poly"]
        rows = poly_table[poly_table["name"].str.startswith("ground_strap_")]
        if len(rows) == 0:
            raise ValueError("No ground_strap_* geometry found.")
        minx, miny, maxx, maxy = rows.geometry.total_bounds
        return minx, miny, maxx, maxy  # already mm

    if component_name not in design.components:
        raise ValueError(f"'{component_name}' not in design.components: "
                          f"{list(design.components.keys())}")
    comp_id = design.components[component_name].id

    bounds_list = []
    for table_name, table in design.qgeometry.tables.items():
        if "component" not in table.columns:
            continue
        rows = table[table["component"] == comp_id]
        if len(rows) == 0:
            continue
        bounds_list.append(rows.geometry.total_bounds)
        print(f"    (found {len(rows)} row(s) for '{component_name}' in "
              f"qgeometry table '{table_name}')")

    if not bounds_list:
        raise ValueError(f"No qgeometry found for component '{component_name}' "
                          f"in ANY table. Tables checked: "
                          f"{list(design.qgeometry.tables.keys())}")

    import numpy as np
    bounds_arr = np.array(bounds_list)  # each row: [minx,miny,maxx,maxy]
    minx, miny = bounds_arr[:, 0].min(), bounds_arr[:, 1].min()
    maxx, maxy = bounds_arr[:, 2].max(), bounds_arr[:, 3].max()
    return minx, miny, maxx, maxy  # mm - matches .pvd file coordinate units


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ground-pos", type=parse_ground_pos, default=5.906)
    p.add_argument("--total-length-mm", type=float, default=None)
    p.add_argument("--component", default=None,
                    help="Component name, e.g. resonator1, TL, LP1, LP2, "
                         "ground_strap. See --list-components.")
    p.add_argument("--list-components", action="store_true")
    args = p.parse_args()

    total_length_mm = args.total_length_mm if args.total_length_mm is not None else TOTAL_LENGTH_MM
    design = build_design(ground_pos_mm=args.ground_pos, total_length_mm=total_length_mm)
    design.rebuild()

    if args.list_components:
        print("Design components:", list(design.components.keys()))
        print("qgeometry tables:", list(design.qgeometry.tables.keys()))
        for table_name, table in design.qgeometry.tables.items():
            print(f"  '{table_name}': {len(table)} rows, "
                  f"columns={list(table.columns)}")
        return

    if not args.component:
        print("ERROR: --component is required (or use --list-components).")
        sys.exit(1)

    xmin, ymin, xmax, ymax = get_component_bbox_mm(design, args.component)
    print(f"\nBounding box for '{args.component}' (meters):")
    print(f"  X: [{xmin:.6e}, {xmax:.6e}]")
    print(f"  Y: [{ymin:.6e}, {ymax:.6e}]")
    print(f"\nPaste into field_viz_by_component.py (pvpython):")
    print(f"  --xmin {xmin:.6e} --ymin {ymin:.6e} --xmax {xmax:.6e} --ymax {ymax:.6e}")


if __name__ == "__main__":
    main()
