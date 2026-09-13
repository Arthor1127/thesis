"""Minimal parser for the subset of gds2palace's XML stackup schema this
pipeline needs: a single substrate material and a single metal layer.

Reimplemented against stdlib xml.etree.ElementTree rather than depending on
gds2palace itself -- gds2palace requires Python 3.12+ (PEP 701 f-strings,
confirmed at util_elmer.py:487 of the installed package) and cannot be
installed alongside SQDMetal, which is pinned to <3.12. The schema fields
used here (Material/Permittivity, Dielectric/Material/Thickness,
Layer/Layer/Type) were confirmed by reading gds2palace's
util_stackup_reader.py directly (see notebooks/design/gds_stackup_spike/
COMPARISON_SQDMETAL_VS_GDS2PALACE.md).

Multi-layer/multi-material stackups are not supported -- same scope limit
as the Phase 1 spike.
"""
from dataclasses import dataclass
import xml.etree.ElementTree as ET


@dataclass
class StackupInfo:
    material_name: str
    permittivity: float
    chip_thickness_m: float
    metal_gds_layer: int


def resolve_stackup(xml_path: str) -> StackupInfo:
    root = ET.parse(xml_path).getroot()

    materials = {m.get("Name"): float(m.get("Permittivity", 1)) for m in root.iter("Material")}

    dielectrics = list(root.iter("Dielectric"))
    assert len(dielectrics) == 1, "this pipeline only supports a single-substrate stackup"
    substrate = dielectrics[0]
    material_name = substrate.get("Material")
    assert material_name in materials, f"Dielectric references unknown Material '{material_name}'"

    layers = list(root.iter("Layer"))
    assert len(layers) == 1, "this pipeline only supports a single metal layer"
    metal_layer = int(layers[0].get("Layer"))

    return StackupInfo(
        material_name=material_name,
        permittivity=materials[material_name],
        chip_thickness_m=float(substrate.get("Thickness")),
        metal_gds_layer=metal_layer,
    )
