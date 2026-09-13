"""Thin shim -- the real code now lives in the ``qmview`` package (installed in
the qcg-quantum-design env). Kept so existing notebooks that do
``from render_geometry import render, sanity_checks`` still work.

    from qmview import render, sanity_checks, preview, ground_pieces
"""
from qmview import (  # noqa: F401
    render, preview, sanity_checks, ground_pieces, design_bounds_mm, poly_patches,
    poly_patches as _poly_patches,
)
