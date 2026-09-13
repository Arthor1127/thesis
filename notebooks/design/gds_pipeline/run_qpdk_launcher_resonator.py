"""Driver for the Phase 3 qpdk launcher+resonator eigenmode run. Reads
qpdk_run_config.json (written by make_qpdk_launcher_resonator.py) and calls
gds_stackup_pipeline.build.run_gds_eigenmode with a single geometry-only
CPW-edge port at the launcher's probe end.

Local mesh-preview use (no simulation, safe on this machine):
    conda run -n qcg-general-design python run_qpdk_launcher_resonator.py --no-run

Cluster use (inside submit_qpdk_launcher_resonator.sh):
    python -u run_qpdk_launcher_resonator.py --num-cpus "$NSLOTS" \\
        --palace-bin "$PALACE_BIN" --run
"""
import argparse
import json
import os
import sys

import matplotlib

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "packages", "gds-stackup-pipeline", "src"))
from gds_stackup_pipeline import run_gds_eigenmode

# Must happen AFTER importing gds_stackup_pipeline, not before: that import
# pulls in SQDMetal's Eigenmode_Simulation -> PVDVTU_Viewer -> `import
# pyvista`, and pyvista resets the matplotlib backend to a Qt one as a side
# effect of initializing VTK -- confirmed directly (instrumented trace:
# backend was "Agg" right after an earlier matplotlib.use() call, flipped to
# "QtAgg" the moment gds_stackup_pipeline was imported, nothing else in
# between). MPLBACKEND=Agg in the job's environment doesn't survive this
# either, since it's pyvista actively switching the backend, not a case of
# the env var being ignored. Forcing it again here, after the import that
# causes the reset, is what actually sticks.
matplotlib.use("Agg", force=True)


def _build_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=HERE)
    p.add_argument("--sim-name", default="qpdk_launcher_resonator")
    p.add_argument("--num-cpus", type=int, default=32)
    p.add_argument("--palace-bin", default="palace")
    p.add_argument("--run", dest="run", action="store_true", default=True)
    p.add_argument("--no-run", dest="run", action="store_false")
    return p


def main():
    args = _build_parser().parse_args()

    with open(os.path.join(HERE, "qpdk_run_config.json")) as f:
        config = json.load(f)

    outdir = args.outdir if args.outdir.endswith("/") else args.outdir + "/"
    os.makedirs(outdir, exist_ok=True)

    port = {
        "type": "cpw_edge",
        "pt_near_centre_end": tuple(config["port_pt_near_centre_end_m"]),
        "len_launch": config["port_len_launch_m"],
    }

    result = run_gds_eigenmode(
        # Resolved relative to this script's own location, not config's
        # baked-in absolute paths -- those were written at generation time
        # on whatever machine ran make_qpdk_launcher_resonator.py, and don't
        # exist once this deploys to a different machine (e.g. the cluster).
        gds_path=os.path.join(HERE, os.path.basename(config["gds_path"])),
        stackup_xml_path=os.path.join(HERE, os.path.basename(config["stackup_xml_path"])),
        sim_name=args.sim_name,
        sim_parent_directory=outdir,
        ports=[port],
        max_feature_size=config["max_feature_size_m"],
        segment_rectangle=tuple(config["segment_rectangle_m"]),
        user_options={"num_cpus": args.num_cpus, "palace_dir": args.palace_bin},
        run=args.run,
    )

    print(result["status"].upper())
    if result["error"]:
        print(result["error"])
        sys.exit(1)


if __name__ == "__main__":
    main()
