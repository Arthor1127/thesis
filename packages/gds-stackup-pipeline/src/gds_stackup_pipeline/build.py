"""Run a SQDMetal PALACE eigenmode simulation from a raw GDS layout + an XML
material stackup, via SQDMetal's native `gds_design=` input branch (proven
working in notebooks/design/gds_stackup_spike/).

Mirrors the call sequence used for QDesign input in
notebooks/design/squid_array_coupling/build_and_eigenmode.py:233-307, minus
the QDesign-only calls (add_ground_plane, fine_mesh_components,
create_port_JosephsonJunction) which don't work with GDS input (see
COMPARISON_SQDMETAL_VS_GDS2PALACE.md).

Ports are geometry-only, decided convention (b): the caller supplies port
geometry directly (no GDS-layer port convention, no sidecar port file).
"""
import json
import os

from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
from SQDMetal.Utilities.Materials import Material

from .stackup import resolve_stackup


def _json_safe(obj):
    """Recursively convert numpy arrays/scalars into plain Python types so
    json.dump doesn't choke on them. Same pattern as build_and_eigenmode.py's
    _json_safe() -- numpy complex .tolist() still isn't JSON serializable,
    so that result must be re-checked."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, complex):
        return {"real": obj.real, "imag": obj.imag}
    if hasattr(obj, "tolist"):
        return _json_safe(obj.tolist())
    return obj


def _add_port(sim, port):
    port = dict(port)
    port_type = port.pop("type")
    if port_type == "waveport":
        plane = port.pop("plane")
        sim.create_waveport_on_boundary(plane, **port)
    elif port_type == "cpw_edge":
        pt = port.pop("pt_near_centre_end")
        len_launch = port.pop("len_launch")
        sim.create_port_CPW_via_edge_point(pt, len_launch, **port)
    else:
        assert False, f"Unknown port type '{port_type}' (expected 'waveport' or 'cpw_edge')"


def run_gds_eigenmode(gds_path, stackup_xml_path, sim_name, sim_parent_directory,
                       ports, fine_mesh_rect=None, max_feature_size=None,
                       segment_rectangle=None, user_options=None, run=True):
    """
    Parameters
    ----------
    gds_path : str
        Path to the GDS layout file.
    stackup_xml_path : str
        Path to the gds2palace-schema XML stackup (single substrate, single
        metal layer -- see stackup.resolve_stackup).
    sim_name : str
        Simulation name, passed through to PALACE_Eigenmode_Simulation.
    sim_parent_directory : str
        Output directory (must end with a slash), passed through.
    ports : list[dict]
        One dict per port. Each needs a 'type' key, 'waveport' or 'cpw_edge':
          - 'waveport': {'type': 'waveport', 'plane': 'z_pos', 'x':, 'y':, 'size_x':, 'size_y':}
            (or the x/y/z + size_* combination create_waveport_on_boundary expects for other planes)
          - 'cpw_edge': {'type': 'cpw_edge', 'pt_near_centre_end': (x, y), 'len_launch': float}
    fine_mesh_rect : tuple[float, float, float, float], optional
        (x1, y1, x2, y2) rectangle to fine-mesh, in metres, via
        fine_mesh_in_rectangle(). This registers a filled surface region --
        GMSH's distance field for a region is zero *everywhere inside* it,
        not just near the actual metal, so the whole rectangle interior gets
        meshed at mesh_min uniformly; taper only happens outside the
        rectangle. Fine for a rectangle deliberately sized to just the
        feature that needs it, wrong (wastefully fine, slow to solve) if
        sized to cover most/all of the chip -- confirmed directly on the
        qpdk launcher+resonator run, where sizing it to the metal's own
        bbox produced a uniformly fine mesh with no taper anywhere inside
        the circuit's own footprint. Prefer max_feature_size for CPW-style
        geometry (see below). Exactly one of fine_mesh_rect/max_feature_size
        is required -- SQDMetal's mesh builder errors with zero fine-mesh
        regions registered, and there's no QDesign to infer a default region
        from for GDS input.
    max_feature_size : float, optional
        Alternative to fine_mesh_rect, via fine_mesh_features(). Geometry-aware:
        auto-detects actual narrow features (a CPW trace/gap) below this
        width, in metres, and fine-meshes only those -- a wide launcher pad
        or ground-plane bulk stays coarse. The right choice for CPW-style
        designs; pick a value comfortably above the trace/gap width and
        comfortably below any wide pad width.
    segment_rectangle : tuple[float, float, float, float], optional
        (x1, y1, x2, y2) substrate/chip rectangle, in metres. Without it,
        GeomGDS infers the chip bounding box exactly from the metal's own
        bbox -- every metal edge then touches the substrate boundary
        exactly, which for some geometries produces a degenerate
        (near-zero-length) polygon edge in SQDMetal's substrate-minus-metal
        gap computation, and GMSH refuses to draw it ("Could not create
        line"). Confirmed with the qpdk launcher+resonator geometry: the
        fix isn't a metal-simplification threshold (large enough thresholds
        to avoid the crash also bridge the real 6um CPW gap shut, visibly
        corrupting the design) -- it's giving the substrate a few
        micrometres of margin beyond the metal's own bbox so no edge
        coincides exactly. Pass this whenever a run hits that error.
    user_options : dict, optional
        Passed through to PALACE_Eigenmode_Simulation's user_options.
    run : bool
        If True (default), calls sim.run() after prepare_simulation() --
        this is the one step with real side effects (submits/executes the
        actual Palace solve). Set False to stop after mesh generation only.

    Returns
    -------
    dict with 'status' ('pass'/'fail'), 'mesh_file' (path or None on
    failure), 'eigenfrequencies' (the real/imag eigenfrequency + Q array
    PALACE_Eigenmode_Simulation.run() returns, or None if run=False or on
    failure), and 'error' (traceback string or None).
    """
    stackup = resolve_stackup(stackup_xml_path)

    try:
        Material(stackup.material_name)
    except AssertionError as e:
        return {"status": "fail", "mesh_file": None, "eigenfrequencies": None,
                "error": f"No SQDMetal material match for '{stackup.material_name}': {e}"}

    result = {"status": "fail", "mesh_file": None, "eigenfrequencies": None, "error": None}
    try:
        sim = PALACE_Eigenmode_Simulation(
            name=sim_name,
            sim_parent_directory=sim_parent_directory,
            mode="simPC",
            meshing="GMSH",
            create_files=True,
            gds_design=gds_path,
            dielectric_material=stackup.material_name,
            gds_chip_thickness=stackup.chip_thickness_m,
            segment_rectangle=segment_rectangle,
            user_options=user_options or {},
        )
        sim.add_metallic(stackup.metal_gds_layer)

        assert (fine_mesh_rect is None) != (max_feature_size is None), (
            "exactly one of fine_mesh_rect or max_feature_size is required -- "
            "SQDMetal's mesh builder errors with zero fine-mesh regions "
            "registered, and there's no QDesign to infer a default region "
            "from for GDS input."
        )
        if max_feature_size is not None:
            sim.fine_mesh_features(max_feature_size)
        else:
            sim.fine_mesh_in_rectangle(*fine_mesh_rect)

        for port in ports:
            _add_port(sim, port)

        sim.prepare_simulation()
        result["mesh_file"] = os.path.join(sim_parent_directory, sim_name, f"{sim_name}.msh")
        if run:
            result["eigenfrequencies"] = _json_safe(sim.run())
        result["status"] = "pass"
    except Exception:
        import traceback
        result["error"] = traceback.format_exc()

    with open(os.path.join(sim_parent_directory, f"{sim_name}_result.json"), "w") as f:
        json.dump(_json_safe(result), f, indent=2)

    return result
