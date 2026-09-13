"""Backwards-compatible shim -- the colormaps moved to :mod:`thesis_cmaps`.

``import qcg_cmaps`` still registers the ``qcg`` / ``qcg_r`` colormaps (now
legacy aliases of ``thesis_cmap`` / ``thesis_cmap_r``).  New code should
``import thesis_cmaps`` and use ``cmap="thesis_cmap"`` /
``cmap="thesis_cmap_diverging"``.
"""

from __future__ import annotations

from thesis_cmaps import (  # noqa: F401
    THESIS_PALETTE,
    register,
    thesis_cmap,
    thesis_cmap_diverging,
    thesis_cmap_diverging_r,
    thesis_cmap_r,
)
from thesis_cmaps import THESIS_PALETTE as QCG_PALETTE  # noqa: F401
from thesis_cmaps import thesis_cmap as qcg_cmap  # noqa: F401
from thesis_cmaps import thesis_cmap_r as qcg_cmap_r  # noqa: F401
