"""Spike: prove a GDS + gds2palace-style XML stackup (resolved to JSON by
parse_stackup_to_json.py) can be adapted into SQDMetal's existing native
GDS ingestion + Gmsh-building code, with no changes to SQDMetal itself.

Run under an env with SQDMetal (e.g. qcg-general-design, cloned from
qcg-quantum-design). Requires resolved_stackup.json to already exist --
produced by running parse_stackup_to_json.py under an env with gds2palace
(e.g. qcg-gds-design) first. The two can't share one env: gds2palace needs
Python 3.12+, SQDMetal needs <3.12.
"""
import json
import os
import shutil
import traceback

import gdstk

from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.Utilities.Materials import Material

HERE = os.path.dirname(os.path.abspath(__file__))
GDS_PATH = os.path.join(HERE, "test.gds")
SIM_DIR = os.path.join(HERE, "spike_sim_output") + "/"


def make_test_gds(path):
    lib = gdstk.Library()  # default unit=1e-6 (coordinates given in microns)
    cell = lib.new_cell("TOP")
    # a single 100x100 micron metal pad on layer 1 -- all the geometry this spike needs
    cell.add(gdstk.rectangle((-50, -50), (50, 50), layer=1))
    lib.write_gds(path)


def resolve_material(material_name):
    """Reuse SQDMetal's own Material lookup rather than hand-copying its table."""
    try:
        mat = Material(material_name)
        return material_name, True
    except AssertionError:
        return None, False


def run_spike():
    with open(os.path.join(HERE, "resolved_stackup.json")) as f:
        stackup = json.load(f)

    make_test_gds(GDS_PATH)

    resolved_name, exact_match = resolve_material(stackup["material_name"])
    assert resolved_name is not None, (
        f"No SQDMetal material match for '{stackup['material_name']}' -- "
        "spike stops here, this is a real adapter gap, not a bug to paper over."
    )

    if os.path.isdir(SIM_DIR):
        shutil.rmtree(SIM_DIR)
    os.makedirs(SIM_DIR, exist_ok=True)

    result = {
        "status": "fail",
        "resolved_material": resolved_name,
        "material_exact_match": exact_match,
        "error": None,
    }

    try:
        sim = PALACE_Eigenmode_Simulation(
            name="gds_stackup_spike",
            sim_parent_directory=SIM_DIR,
            mode="simPC",
            meshing="GMSH",
            create_files=True,
            gds_design=GDS_PATH,
            dielectric_material=resolved_name,
            gds_chip_thickness=stackup["chip_thickness_m"],
        )
        sim.add_metallic(stackup["metal_gds_layer"])
        # GMSH_Mesh_Builder.build_mesh() requires at least one fine-mesh field;
        # fine_mesh_components() is QDesign-only, so use the geometry-only
        # fine_mesh_in_rectangle() instead -- covers the metal pad footprint.
        sim.fine_mesh_in_rectangle(-50e-6, -50e-6, 50e-6, 50e-6)
        sim.create_waveport_on_boundary(
            "z_pos", x=0, y=0, size_x=20e-6, size_y=20e-6
        )
        sim.prepare_simulation()
        result["status"] = "pass"
    except Exception:
        result["error"] = traceback.format_exc()

    with open(os.path.join(HERE, "spike_result.json"), "w") as f:
        json.dump(result, f, indent=2)

    print(result["status"].upper())
    if result["error"]:
        print(result["error"])
    return result


if __name__ == "__main__":
    run_spike()
