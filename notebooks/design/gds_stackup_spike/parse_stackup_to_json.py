"""Parse a gds2palace-schema XML stackup into the small JSON handoff shape
SQDMetal's side of the spike needs. Run under an env that has gds2palace
(e.g. qcg-gds-design) -- gds2palace requires Python 3.12+ (PEP 701 f-strings
in util_elmer.py), SQDMetal requires <3.12, so the two can't share one env.
"""
import json
import sys

import gds2palace.util_stackup_reader as usr


def resolve_stackup(xml_path):
    materials, dielectrics, metals = usr.read_substrate(xml_path)

    assert len(dielectrics.dielectrics) == 1, "spike only handles a single-substrate stackup"
    assert len(metals.metals) == 1, "spike only handles a single metal layer"

    substrate = dielectrics.dielectrics[0]
    material = materials.get_by_name(substrate.material)
    metal = metals.metals[0]

    return {
        "material_name": material.name,
        "permittivity": material.eps,
        "chip_thickness_m": substrate.thickness,
        "metal_gds_layer": int(metal.layernum),
    }


def demo():
    result = resolve_stackup("test_stackup.xml")
    assert result["material_name"] == "Silicon"
    assert abs(result["permittivity"] - 11.45) < 1e-9
    assert abs(result["chip_thickness_m"] - 500e-6) < 1e-12
    assert result["metal_gds_layer"] == 1
    print("demo OK:", result)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--demo":
        demo()
        sys.exit(0)

    resolved = resolve_stackup("test_stackup.xml")
    with open("resolved_stackup.json", "w") as f:
        json.dump(resolved, f, indent=2)
    print("wrote resolved_stackup.json:", resolved)
