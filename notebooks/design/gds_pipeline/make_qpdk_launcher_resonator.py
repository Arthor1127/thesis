"""Generate a real qpdk geometry (launcher + coupled quarter-wave CPW
resonator) and its companion stackup XML + run config, for the Phase 3
cluster eigenmode run. Geometry generation only -- no simulation, safe to
run locally.

Run in the qcg-gds-design conda env (has gdsfactory + qpdk):
    conda run -n qcg-gds-design python make_qpdk_launcher_resonator.py

Outputs (next to this script):
    qpdk_launcher_resonator.gds  -- the layout (metal on GDS layer (1, 0))
    qpdk_stackup.xml             -- material stackup, values from qpdk.tech
    qpdk_run_config.json         -- port coordinate + fine-mesh rectangle,
                                     read by run_qpdk_launcher_resonator.py
"""
import json
import math
import os

import gdsfactory as gf
from qpdk import PDK
from qpdk.cells.launcher import launcher
from qpdk.cells.resonator import resonator_coupled
from qpdk.tech import LAYER, material_properties

HERE = os.path.dirname(os.path.abspath(__file__))
UM = 1e-6

SUBSTRATE_THICKNESS_UM = 500  # qpdk.tech.get_layer_stack()'s "Substrate" LayerLevel
FINE_MESH_MARGIN_UM = 50
SUBSTRATE_MARGIN_UM = 200  # matches qpdk.cells.chip.chip_edge()'s own default frame width


def build_component():
    """qpdk's coplanar_waveguide cross-section only draws the center trace on
    M1_DRAW plus two etch-gap marker strips on M1_ETCH (confirmed by reading
    qpdk.tech.coplanar_waveguide's cross-section sections directly) -- it
    does NOT include a ground-plane fill, since qpdk expects that to come
    from a separate full-chip layer added elsewhere. A launcher+resonator
    component built from just these cells is therefore an ungrounded
    floating trace, not a real CPW: SQDMetal's create_CPW_feed_via_point
    (which looks for a center conductor flanked by two ground pads) fails
    on it with "CPW couldn't be formed" since there's no ground metal
    anywhere nearby to find.

    Fix applied here: draw an explicit ground-plane rectangle over the
    assembly's bounding box (plus margin) and boolean-subtract the M1_ETCH
    gap layer from it, then union with the M1_DRAW trace -- a real grounded
    CPW system SQDMetal can actually port and mesh, and a closer match to
    how these designs are used in practice (ground fill is not optional).

    The west (launcher/probe) edge is deliberately NOT padded with margin --
    it stays exactly at the launcher's own outer pad edge (x=0), both because
    a real launcher sits flush with the chip edge for probe/wirebond access
    (solid ground metal in front of it would be wrong), and because the
    launcher's "waveport" reference point sits exactly on that edge: padding
    it would bury the port point inside solid ground metal instead of
    leaving it on a boundary, which silently breaks
    create_port_CPW_via_edge_point's edge-finding (confirmed directly: with
    padding, poly_all.contains(Point(0,0)) is True, and the port constructor
    finds the ground rectangle's own flat outer edge -- with no CPW corner
    structure nearby -- instead of the launcher's actual pad edge).
    """
    PDK.activate()
    raw = gf.Component()
    launcher_ref = raw.add_ref(launcher())
    resonator_ref = raw.add_ref(resonator_coupled())
    resonator_ref.connect("coupling_o1", launcher_ref.ports["o1"])

    margin = FINE_MESH_MARGIN_UM
    ground = gf.Component()
    ground.add_polygon(
        [(raw.xmin, raw.ymin - margin), (raw.xmax + margin, raw.ymin - margin),
         (raw.xmax + margin, raw.ymax + margin), (raw.xmin, raw.ymax + margin)],
        layer=LAYER.M1_DRAW,
    )
    ground_minus_gap = gf.boolean(ground, raw, "not", layer=LAYER.M1_DRAW,
                                   layer1=LAYER.M1_DRAW, layer2=LAYER.M1_ETCH)
    grounded = gf.boolean(ground_minus_gap, raw, "or", layer=LAYER.M1_DRAW,
                           layer1=LAYER.M1_DRAW, layer2=LAYER.M1_DRAW)
    return grounded, launcher_ref, raw


def write_stackup_xml(path):
    si = material_properties["Si"]
    with open(path, "w") as f:
        f.write(f"""<?xml version="1.0" encoding="UTF-8"?>
<Stackup>
    <Materials>
        <Material Name="Silicon" Type="Dielectric" Permittivity="{si["relative_permittivity"]}" DielectricLossTangent="{si["loss_tangent"]}" />
    </Materials>
    <Dielectrics>
        <Dielectric Name="Substrate" Material="Silicon" Zmin="{-SUBSTRATE_THICKNESS_UM * UM}" Zmax="0" Thickness="{SUBSTRATE_THICKNESS_UM * UM}" />
    </Dielectrics>
    <Layers>
        <Layer Name="M1" Layer="{LAYER.M1_DRAW[0]}" Type="CONDUCTOR" Material="Silicon" Zmin="0" Zmax="0" />
    </Layers>
</Stackup>
""")


def main():
    c, launcher_ref, raw = build_component()

    gds_path = os.path.join(HERE, "qpdk_launcher_resonator.gds")
    # with_metadata=False: gdsfactory's default write_gds() prepends a
    # "$$$CONTEXT_INFO$$$" metadata pseudo-cell as cells[0], pushing the
    # real top cell to index 1 -- SQDMetal's GeomGDS defaults to
    # gds_cell_index=0 (matches plain gdstk/gdspy-written files, which don't
    # have this cell), so it would silently read zero real polygons.
    c.write_gds(gds_path, with_metadata=False)
    print(f"wrote {gds_path} -- bbox ({c.xmin:.1f}, {c.ymin:.1f}) to "
          f"({c.xmax:.1f}, {c.ymax:.1f}) um, layers {sorted(c.layers)}")

    xml_path = os.path.join(HERE, "qpdk_stackup.xml")
    write_stackup_xml(xml_path)
    print(f"wrote {xml_path}")

    # create_port_CPW_via_edge_point needs its reference point strictly off
    # the metal boundary (inside or outside), not exactly on it: the
    # algorithm computes a perpendicular direction as
    # (nearest_boundary_point - reference_point), normalized -- if the
    # reference point already sits exactly on the boundary, that vector is
    # zero and normalizing it produces NaN, which then makes every
    # downstream "is this corner collinear" check silently false (0 corners
    # found on either side, not any assertion until the final count check).
    # Confirmed this is exactly what happens with qpdk's "waveport" port
    # coordinate used raw: it sits exactly on the launcher pad's edge, and
    # whether the resulting exact-zero-distance case "accidentally" works
    # depends on which GEOS version computes shortest_line() with enough
    # floating-point noise to avoid a true zero -- worked locally (GEOS
    # 3.10.6), failed on the cluster (GEOS 3.13.1). Not a GEOS-version bug
    # to work around: the input was never a valid "near the edge, not on
    # it" reference point. Fixed by insetting a few um into the metal,
    # opposite the port's own outward-facing orientation.
    port_xy_um = launcher_ref.ports["waveport"].center
    inset_um = 5.0
    orientation_rad = math.radians(launcher_ref.ports["waveport"].orientation)
    port_xy_um = (
        port_xy_um[0] - inset_um * math.cos(orientation_rad),
        port_xy_um[1] - inset_um * math.sin(orientation_rad),
    )
    config = {
        "gds_path": gds_path,
        "stackup_xml_path": xml_path,
        "port_pt_near_centre_end_m": [port_xy_um[0] * UM, port_xy_um[1] * UM],
        "port_len_launch_m": 20e-6,
        # fine_mesh_features (geometry-aware: auto-detects narrow features
        # below this width) instead of a fine_mesh_rect covering the whole
        # circuit -- a filled-rectangle region meshes uniformly fine
        # *everywhere inside it* (GMSH's distance field for a region is
        # zero throughout the interior, not just near the metal), which
        # wastefully fine-meshes the whole layout including the wide
        # launcher pad and ground plane, confirmed directly on this design.
        # qpdk's CPW cross-section is 10um trace / 6um gap (tech.py,
        # coplanar_waveguide defaults); the launcher's big pad is 200um
        # wide -- 20um sits comfortably above the former, well below the
        # latter.
        "max_feature_size_m": 20e-6,
        # SUBSTRATE_MARGIN_UM of substrate beyond the metal's own bbox on
        # every side. Two independent reasons this needs to be a real
        # physical margin, not a token amount: (1) with zero/near-zero
        # margin every metal edge touches the substrate boundary exactly,
        # producing a degenerate near-zero-length polygon edge in SQDMetal's
        # substrate-minus-metal gap computation that GMSH refuses to draw
        # ("Could not create line") -- confirmed a metal-simplification
        # threshold large enough to dodge that crash is also large enough
        # to bridge the real 6um CPW gap shut, visibly corrupting the
        # meander; (2) simulation boundary proximity is a real physical
        # effect -- the airbox boundary condition perturbs the eigenmode
        # fields if it sits too close to the circuit, so this margin also
        # needs to be generous on its own physical merits, not just enough
        # to avoid the GMSH crash. This
        # substrate-margin fix avoids the crash without touching the gap.
        "segment_rectangle_m": [
            (c.xmin - SUBSTRATE_MARGIN_UM) * UM,
            (c.ymin - SUBSTRATE_MARGIN_UM) * UM,
            (c.xmax + SUBSTRATE_MARGIN_UM) * UM,
            (c.ymax + SUBSTRATE_MARGIN_UM) * UM,
        ],
    }
    config_path = os.path.join(HERE, "qpdk_run_config.json")
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)
    print(f"wrote {config_path}: {config}")


if __name__ == "__main__":
    main()
