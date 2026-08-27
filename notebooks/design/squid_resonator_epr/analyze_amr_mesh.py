"""
Parse an MFEM mesh (v1.0) file written by Palace's AMR and bin element
centroids into named geometric regions, to see WHERE AMR concentrated
refinement without needing to run another (expensive) solve.

Regions reuse the pad/squid-arm bounding boxes derived in the
2026-08-24/25 chat from squid_options geometry (pad_radius=90um,
theta_2=180, delta_angle=7.5, l1=l2=300um, w1=6um, g1=4um), plus a
readout-region box and an "other" bucket for everything else. Coordinates
in the .mesh file are in the model's L0 units (mm, per L0=0.001 in the
Palace json), so region boxes are converted from meters to mm here.
"""
import sys
import numpy as np

MESH_PATH = sys.argv[1]

PAD_RADIUS_MM = 0.090
PAD_MARGIN_MM = 0.020
SQUID_ARM_BOX_MM = (-0.710, -0.025, -0.080, 0.025)   # x1,y1,x2,y2

# Approximate readout resonator bounding box - it runs from the qubit's
# north side out across the chip; use a generous box that starts past
# the qubit pad going away from the SQUID arms (+y-ish region) and
# excludes the qubit-local area. This is approximate/coarse since we
# only need "is this element near the resonator trace" for a first look.
READOUT_BOX_MM = (-1.2, 0.15, 1.2, 2.0)


def parse_mfem_mesh(path):
    with open(path) as f:
        line = f.readline()
        while not line.startswith("elements"):
            line = f.readline()
        n_elem = int(f.readline())
        elem_verts = np.empty((n_elem, 4), dtype=np.int64)
        for i in range(n_elem):
            parts = f.readline().split()
            # format: <attr> <geom_type> <v0> <v1> <v2> <v3>
            elem_verts[i] = [int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5])]

        line = f.readline()
        while not line.startswith("vertices"):
            line = f.readline()
        n_vert = int(f.readline())
        # Curved/high-order mesh: vertex coords come from a "nodes"
        # GridFunction (H1_3D_P1, VDim=3, Ordering: 1 = interleaved
        # x,y,z per line), not a flat coordinate list directly after the
        # vertex count - see 2026-08-25 chat. Skip the preamble up to
        # and including the blank line after "Ordering: <n>".
        line = f.readline()
        while not line.startswith("Ordering:"):
            line = f.readline()
        f.readline()  # blank line separating the GridFunction header from data
        verts = np.empty((n_vert, 3), dtype=np.float64)
        for i in range(n_vert):
            verts[i] = [float(x) for x in f.readline().split()]

    return elem_verts, verts


def main():
    print(f"Parsing {MESH_PATH} ...", flush=True)
    elem_verts, verts = parse_mfem_mesh(MESH_PATH)
    print(f"  {len(elem_verts)} elements, {len(verts)} vertices", flush=True)

    centroids = verts[elem_verts].mean(axis=1)  # (n_elem, 3)
    cx, cy = centroids[:, 0], centroids[:, 1]

    counts = {"pad": 0, "squid_arm": 0, "readout": 0, "other": 0}
    # Vectorized region binning
    x1, y1, x2, y2 = SQUID_ARM_BOX_MM
    in_squid = (cx >= x1) & (cx <= x2) & (cy >= y1) & (cy <= y2)
    in_pad = (~in_squid) & (cx * cx + cy * cy <= (PAD_RADIUS_MM + PAD_MARGIN_MM) ** 2)
    rx1, ry1, rx2, ry2 = READOUT_BOX_MM
    in_readout = (~in_squid) & (~in_pad) & (cx >= rx1) & (cx <= rx2) & (cy >= ry1) & (cy <= ry2)
    in_other = ~(in_squid | in_pad | in_readout)

    counts["squid_arm"] = int(in_squid.sum())
    counts["pad"] = int(in_pad.sum())
    counts["readout"] = int(in_readout.sum())
    counts["other"] = int(in_other.sum())

    total = len(elem_verts)
    print("\nElement distribution by region:")
    for region, n in counts.items():
        print(f"  {region:>10}: {n:>10}  ({100*n/total:.2f}%)")
    print(f"  {'TOTAL':>10}: {total:>10}")


if __name__ == "__main__":
    main()
