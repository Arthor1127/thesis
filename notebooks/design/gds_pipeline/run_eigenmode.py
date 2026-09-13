"""CLI driver: run a SQDMetal PALACE eigenmode simulation from a raw GDS
layout + an XML material stackup. Thin wrapper around
gds_stackup_pipeline.build.run_gds_eigenmode -- see that function's
docstring and packages/gds-stackup-pipeline/README.md for details.

There's no qiskit_metal "build" step here (the GDS is externally authored),
so unlike build_and_eigenmode.py there's no build_design.py sibling -- every
knob is a CLI flag.

Example (Phase 1 spike fixtures). Note --fine-mesh-rect needs plain decimal
form, not scientific notation -- argparse's negative-number detection
doesn't recognize "-50e-6" and treats it as an unknown flag. Use
--max-feature-size instead for CPW-style geometry (see its help text):
    conda run -n qcg-general-design python run_eigenmode.py \\
        --gds ../gds_stackup_spike/test.gds \\
        --stackup-xml ../gds_stackup_spike/test_stackup.xml \\
        --outdir /tmp/gds_pipeline_smoke/ \\
        --fine-mesh-rect -0.00005 -0.00005 0.00005 0.00005 --no-run
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "..", "..", "..", "packages", "gds-stackup-pipeline", "src"))
from gds_stackup_pipeline import run_gds_eigenmode


def _build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gds", required=True, help="Path to the GDS layout file.")
    p.add_argument("--stackup-xml", required=True, help="Path to the XML material stackup.")
    p.add_argument("--outdir", required=True, help="Output directory (sim_parent_directory).")
    p.add_argument("--sim-name", default="gds_eigenmode")

    p.add_argument("--fine-mesh-rect", type=float, nargs=4, metavar=("X1", "Y1", "X2", "Y2"),
                    help="Fine-mesh a filled rectangular region (fine_mesh_in_rectangle) -- "
                         "GMSH's distance field for a filled region is zero *everywhere "
                         "inside it*, not just near the metal, so sizing this to cover most "
                         "of the chip meshes the whole interior at mesh_min uniformly with no "
                         "taper (confirmed directly: a launcher+CPW-resonator run with this "
                         "rectangle sized to the metal's own bbox produced a 3x larger mesh "
                         "than the same design meshed with --max-feature-size, and the "
                         "oversized run failed to converge inside a 1hr cluster job while the "
                         "corrected one finished in ~44min). Prefer --max-feature-size for "
                         "CPW-style geometry. Exactly one of the two is required.")
    p.add_argument("--max-feature-size", type=float,
                    help="Fine-mesh only actual narrow features (fine_mesh_features) below "
                         "this width, in metres -- geometry-aware, auto-detects a CPW "
                         "trace/gap and leaves a wide pad or ground-plane bulk coarse. The "
                         "right choice for CPW-style designs; pick a value comfortably above "
                         "the trace/gap width and comfortably below any wide pad width. "
                         "Exactly one of --fine-mesh-rect/--max-feature-size is required.")

    p.add_argument("--port-plane", default="z_pos",
                    choices=["x_pos", "x_neg", "y_pos", "y_neg", "z_pos", "z_neg"])
    p.add_argument("--port-x", type=float, default=0.0)
    p.add_argument("--port-y", type=float, default=0.0)
    p.add_argument("--port-size-x", type=float, default=20e-6)
    p.add_argument("--port-size-y", type=float, default=20e-6)

    p.add_argument("--run", dest="run", action="store_true", default=True,
                    help="Actually run the Palace solve after mesh generation (default).")
    p.add_argument("--no-run", dest="run", action="store_false",
                    help="Stop after mesh generation -- no Palace solve.")
    return p


def main():
    args = _build_parser().parse_args()

    if (args.fine_mesh_rect is None) == (args.max_feature_size is None):
        _build_parser().error("exactly one of --fine-mesh-rect/--max-feature-size is required")

    outdir = args.outdir if args.outdir.endswith("/") else args.outdir + "/"
    os.makedirs(outdir, exist_ok=True)

    # --port-x/--port-y/--port-size-x/--port-size-y are generic "first/second
    # coordinate" and "first/second size" -- which axis they map onto depends
    # on --port-plane, per create_waveport_on_boundary's own convention
    # (the port is coplanar with the chosen air-box boundary).
    port = {"type": "waveport", "plane": args.port_plane}
    if args.port_plane in ("z_pos", "z_neg"):
        port.update(x=args.port_x, y=args.port_y, size_x=args.port_size_x, size_y=args.port_size_y)
    elif args.port_plane in ("y_pos", "y_neg"):
        port.update(x=args.port_x, z=args.port_y, size_x=args.port_size_x, size_z=args.port_size_y)
    else:
        port.update(y=args.port_x, z=args.port_y, size_y=args.port_size_x, size_z=args.port_size_y)

    result = run_gds_eigenmode(
        gds_path=args.gds,
        stackup_xml_path=args.stackup_xml,
        sim_name=args.sim_name,
        sim_parent_directory=outdir,
        ports=[port],
        fine_mesh_rect=tuple(args.fine_mesh_rect) if args.fine_mesh_rect else None,
        max_feature_size=args.max_feature_size,
        run=args.run,
    )

    print(result["status"].upper())
    if result["error"]:
        print(result["error"])
        sys.exit(1)


if __name__ == "__main__":
    main()
