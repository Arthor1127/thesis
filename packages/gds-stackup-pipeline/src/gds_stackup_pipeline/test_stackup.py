"""Assert-based self-check, no framework. Run directly:
    conda run -n qcg-general-design python -m gds_stackup_pipeline.test_stackup

Reuses the Phase 1 spike's fixtures (test.gds, test_stackup.xml) from
notebooks/design/gds_stackup_spike/.
"""
import os

from .stackup import resolve_stackup
from .build import run_gds_eigenmode

FIXTURE_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..",
    "notebooks", "design", "gds_stackup_spike",
)
FIXTURE_DIR = os.path.abspath(FIXTURE_DIR)


def test_resolve_stackup():
    info = resolve_stackup(os.path.join(FIXTURE_DIR, "test_stackup.xml"))
    assert info.material_name == "Silicon"
    assert abs(info.permittivity - 11.45) < 1e-9
    assert abs(info.chip_thickness_m - 500e-6) < 1e-12
    assert info.metal_gds_layer == 1
    print("test_resolve_stackup: PASS")


def test_run_gds_eigenmode_reaches_gmsh(tmp_dir):
    result = run_gds_eigenmode(
        gds_path=os.path.join(FIXTURE_DIR, "test.gds"),
        stackup_xml_path=os.path.join(FIXTURE_DIR, "test_stackup.xml"),
        sim_name="pkg_smoke_test",
        sim_parent_directory=tmp_dir,
        ports=[{"type": "waveport", "plane": "z_pos", "x": 0, "y": 0,
                "size_x": 20e-6, "size_y": 20e-6}],
        fine_mesh_rect=(-50e-6, -50e-6, 50e-6, 50e-6),
        run=False,
    )
    assert result["status"] == "pass", result["error"]
    assert os.path.isfile(result["mesh_file"]), f"expected mesh file at {result['mesh_file']}"
    print("test_run_gds_eigenmode_reaches_gmsh: PASS")


if __name__ == "__main__":
    import shutil
    import tempfile

    test_resolve_stackup()
    tmp_dir = tempfile.mkdtemp(prefix="gds_stackup_pipeline_test_") + "/"
    try:
        test_run_gds_eigenmode_reaches_gmsh(tmp_dir)
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)
    print("ALL PASS")
