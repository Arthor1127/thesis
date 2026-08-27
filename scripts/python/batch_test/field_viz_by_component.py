"""
Clip an existing Palace .pvd field output (eigenmode or driven) down to a
single design component's spatial bounding box, so you can look at (e.g.)
just the resonator's E-field magnitude without the ports/TL/airbox in the
way. Works on the REAL simulation output - no remeshing, no changes to
build_and_eigenmode.py/build_and_driven.py needed.

This is a TWO-STEP workflow because pvpython (ParaView's own bundled
Python) does NOT have qiskit_metal/SQDMetal installed, and your conda
env's Python does not have the paraview module - they can't be combined
in one script run under one interpreter.

IMPORTANT: coordinates throughout are in MILLIMETERS, not meters. Even
though Palace's own console log prints mesh bounds labeled 'm' (meters),
the actual .pvd/.vtu file it writes stores the same raw numbers GMSH/
qiskit-metal used to build the mesh - i.e. millimeters. Confirmed by
comparing ParaView's GetBounds() against the Palace log for the same run
(log said Xmin=-2.880e-03 m; ParaView's real file showed X=[-2.88,2.88]).

Step 1 (regular conda env, has qiskit_metal) - get the bounding box:
    conda activate qcg-quantum-design
    python get_component_bbox.py --ground-pos 5.906 --component resonator1
    # prints e.g.: --xmin -1.234e+00 --ymin ... --xmax ... --ymax ...  (mm)

Step 2 (pvpython, has paraview) - clip and view/screenshot using those
numbers directly, no qiskit_metal import needed:
    pvpython field_viz_by_component.py \
        --pvd ~/sweep_v2_test/eig_5.906mm/outputFiles/paraview/.../electric_field.pvd \
        --xmin=-1.234e+00 --ymin=... --xmax=... --ymax=... \
        --z-center-um 0 --z-margin-um 50 \
        --field-array U_e --output-image resonator1_field.png

Run with --list-arrays (pvpython, just needs --pvd) if unsure of the
exact field array name Palace wrote.
"""
import argparse
import os
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pvd", help="Path to Palace's .pvd output file.")
    p.add_argument("--xmin", type=float, default=None,
                    help="mm - matches the .pvd file's raw coordinate units, "
                         "same as get_component_bbox.py's output. NOT meters.")
    p.add_argument("--ymin", type=float, default=None, help="mm")
    p.add_argument("--xmax", type=float, default=None, help="mm")
    p.add_argument("--ymax", type=float, default=None, help="mm")
    p.add_argument("--field-array", default="U_e",
                    help="Name of the field array to color by (E-field energy "
                         "density is commonly 'U_e' in Palace output; check "
                         "--list-arrays if this doesn't match).")
    p.add_argument("--z-margin-um", type=float, default=50.0,
                    help="Padding (um) added above/below the metal layer in Z. "
                         "Converted to mm internally to match the .pvd file's "
                         "coordinate units (NOT meters - see xmin/ymin help).")
    p.add_argument("--z-center-um", type=float, default=0.0,
                    help="Z center (um) of the metal layer, matching your "
                         "design's substrate/airbox setup. Check the Palace log's "
                         "'Mesh bounding box' Zmin/Zmax if unsure (note: the log "
                         "reports meters, but this flag is um - convert).")
    p.add_argument("--log-scale", action="store_true",
                    help="Use a log color scale instead of linear. Strongly "
                         "recommended for E-field/energy-density quantities - "
                         "these concentrate sharply at conductor edges, so on a "
                         "linear scale nearly the whole volume renders as the "
                         "lowest color (looks blank/flat) while the real "
                         "structure sits in a thin bright sliver you can't see.")
    p.add_argument("--log-scale-min", type=float, default=None,
                    help="Explicit lower bound for the log color scale. Log "
                         "scale can't represent an exact 0, and auto-rescaling "
                         "to the data's true min (often 0.0 for these fields) "
                         "can clamp/behave badly. Run once without this first - "
                         "the console prints a suggested value from the actual "
                         "data's 1st percentile of nonzero values.")
    p.add_argument("--image-resolution", type=int, nargs=2, default=[1200, 1200],
                    help="Screenshot width/height in pixels (default 1200 1200, "
                         "higher than ParaView's small default so thin bright "
                         "field regions near edges are actually visible).")
    p.add_argument("--output-image", default=None,
                    help="If given, saves a screenshot here instead of just "
                         "printing the clip box (still needs pvpython either way).")
    p.add_argument("--list-arrays", action="store_true",
                    help="Print available field arrays in --pvd and exit.")
    args = p.parse_args()

    from paraview.simple import (
        OpenDataFile, Clip, Show, ColorBy, GetActiveViewOrCreate,
        ResetCamera, SaveScreenshot,
    )

    if not args.pvd:
        print("ERROR: --pvd is required.")
        sys.exit(1)

    reader = OpenDataFile(os.path.expanduser(args.pvd))
    reader.UpdatePipeline()

    full_bounds = reader.GetDataInformation().GetBounds()
    print(f"Full dataset bounds (meters): "
          f"X=[{full_bounds[0]:.6e},{full_bounds[1]:.6e}] "
          f"Y=[{full_bounds[2]:.6e},{full_bounds[3]:.6e}] "
          f"Z=[{full_bounds[4]:.6e},{full_bounds[5]:.6e}]")

    if args.list_arrays:
        info = reader.GetDataInformation()
        print("Point arrays:")
        for i in range(info.GetPointDataInformation().GetNumberOfArrays()):
            print("   ", info.GetPointDataInformation().GetArrayInformation(i).GetName())
        print("Cell arrays:")
        for i in range(info.GetCellDataInformation().GetNumberOfArrays()):
            print("   ", info.GetCellDataInformation().GetArrayInformation(i).GetName())
        return

    if None in (args.xmin, args.ymin, args.xmax, args.ymax):
        print("ERROR: --xmin/--ymin/--xmax/--ymax are required (run "
              "get_component_bbox.py under your conda env first to get these).")
        sys.exit(1)

    xmin, ymin, xmax, ymax = args.xmin, args.ymin, args.xmax, args.ymax
    # .pvd file coordinates are in mm (see --xmin help / discovery in
    # get_component_bbox.py), so convert the um Z args to mm too, NOT meters.
    z_center_mm = args.z_center_um * 1e-3
    z_margin_mm = args.z_margin_um * 1e-3
    zmin, zmax = z_center_mm - z_margin_mm, z_center_mm + z_margin_mm

    print(f"Clip box (mm):")
    print(f"   X: [{xmin:.6e}, {xmax:.6e}]")
    print(f"   Y: [{ymin:.6e}, {ymax:.6e}]")
    print(f"   Z: [{zmin:.6e}, {zmax:.6e}]  (from --z-center-um/--z-margin-um - "
          f"widen if the field looks clipped off)")

    clip = Clip(Input=reader)
    clip.ClipType = "Box"
    clip.ClipType.Position = [xmin, ymin, zmin]
    clip.ClipType.Length = [xmax - xmin, ymax - ymin, zmax - zmin]
    clip.Invert = 0  # keep INSIDE the box; flip to 1 if this shows the outside instead
    clip.UpdatePipeline()

    n_points = clip.GetDataInformation().GetNumberOfPoints()
    n_cells = clip.GetDataInformation().GetNumberOfCells()
    print(f"Clip result: {n_points} points, {n_cells} cells.")
    if n_points == 0:
        print("WARNING: clip is EMPTY - nothing survived the box filter. Likely "
              "causes: (1) --xmin/ymin/xmax/ymax don't actually overlap the "
              "mesh's real coordinate range - open the .pvd in the ParaView GUI "
              "and check the Information tab's 'Bounds' for the whole dataset "
              "and compare against the box printed above; (2) --z-center-um/"
              "--z-margin-um miss the metal layer's real Z position - check the "
              "Palace log's 'Mesh bounding box' Zmin/Zmax for this run; "
              "(3) --invert needs flipping.")

    view = GetActiveViewOrCreate("RenderView")

    from paraview.simple import Hide
    Hide(reader, view)  # the reader's own full-dataset display is often left
    # visible by default after OpenDataFile() - since the full dataset spans
    # the whole airbox (much bigger than any single-component clip), it can
    # sit on top of / around the clip and visually swamp it entirely with a
    # single solid color, which is indistinguishable from "clip isn't
    # colored" without explicitly checking for this.

    disp = Show(clip, view)
    disp.SetRepresentationType("Surface")
    ColorBy(disp, ("POINTS", args.field_array))
    disp.RescaleTransferFunctionToDataRange(True)
    ResetCamera(view)

    data_range = clip.PointData[args.field_array].GetRange()
    print(f"'{args.field_array}' data range in this clip: {data_range}")

    # Histogram/percentiles - tells us exactly how concentrated the field is,
    # and gives a real number for a log-scale floor instead of guessing.
    from paraview import servermanager as sm
    from paraview.vtk.numpy_interface import dataset_adapter as dsa
    fetched = sm.Fetch(clip)
    wrapped = dsa.WrapDataObject(fetched)
    values = None
    try:
        import numpy as np
        values = np.asarray(wrapped.PointData[args.field_array]).ravel()
    except Exception as e:
        print(f"    (could not fetch raw values for histogram: {e})")

    suggested_floor = None
    if values is not None and len(values) > 0:
        import numpy as np
        percentiles = [0, 1, 5, 25, 50, 75, 90, 95, 99, 99.9, 100]
        pvals = np.percentile(values, percentiles)
        print(f"    Percentiles of '{args.field_array}' over {len(values)} points:")
        for pct, val in zip(percentiles, pvals):
            print(f"      p{pct:>5}: {val:.6e}")
        nonzero = values[values > 0]
        print(f"    Fraction of points exactly 0: "
              f"{(len(values)-len(nonzero))/len(values)*100:.1f}%")
        if len(nonzero) > 0:
            suggested_floor = float(np.percentile(nonzero, 1))
            print(f"    Suggested log-scale floor (p1 of nonzero values): "
                  f"{suggested_floor:.6e} - pass with --log-scale-min "
                  f"{suggested_floor:.6e} if the plain --log-scale still looks "
                  f"flat (auto-rescale likely used the true min of 0, which "
                  f"log scale can't represent and may silently clamp badly).")

    from paraview.simple import GetColorTransferFunction, Render
    ctf = GetColorTransferFunction(args.field_array)
    if args.log_scale:
        ctf.UseLogScale = 1
        log_min = args.log_scale_min if args.log_scale_min is not None else suggested_floor
        if log_min is not None:
            ctf.RescaleTransferFunction(log_min, data_range[1])
            print(f"    Manually set log-scale range: "
                  f"[{log_min:.6e}, {data_range[1]:.6e}]"
                  f"{' (auto, from p1 of nonzero)' if args.log_scale_min is None else ''}")

    if data_range[0] == data_range[1]:
        print("WARNING: data range is a single constant value - there is no "
              "variation to color by in this region. Either the field is "
              "genuinely uniform/zero here (check you're not clipping into a "
              "region with no field, e.g. entirely inside a PEC conductor), or "
              "something upstream is wrong.")
    else:
        print(f"    (using {'LOG' if args.log_scale else 'LINEAR'} color scale - "
              f"if the image still looks flat/blank on linear, try --log-scale: "
              f"field quantities like this concentrate sharply near conductor "
              f"edges, so most of the clipped volume can legitimately sit near "
              f"the low end of a linear scale.)")

    Render(view)

    print(f"Display color array set to: {disp.ColorArrayName} "
          f"(if this shows ('', '') or doesn't match --field-array, ColorBy "
          f"silently failed - check the exact array name against --list-arrays "
          f"output, and that it's a POINT array not a CELL array).")

    if args.output_image:
        out_path = os.path.expanduser(args.output_image)
        SaveScreenshot(out_path, view, ImageResolution=args.image_resolution)
        print(f"Screenshot saved to: {out_path}")
    else:
        print("No --output-image given - clip/coloring is set up in the active "
              "view; if running pvpython interactively you can inspect it now, "
              "otherwise pass --output-image to save a PNG non-interactively.")


if __name__ == "__main__":
    main()

