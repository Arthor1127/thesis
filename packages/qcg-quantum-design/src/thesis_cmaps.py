"""Thesis colormaps for matplotlib.

Importing this module registers the following colormaps (each with a ``_r``
reverse) with matplotlib, so ``cmap="thesis_cmap"`` works in any plotting call:

``thesis_cmap`` -- sequential, low -> high

    #02566A  dark blue    (absolute low)
    #03AED2  cyan-blue
    #F8DE22  yellow
    #F45B26  orange
    #D12052  magenta-red
    #9E0E3E  deep crimson (absolute high)

``thesis_cmap_diverging`` -- diverging, negative <- 0 -> positive

    #03AED2  cyan-blue   (negative extreme)
    #FFFFFF  white       (zero / midpoint)
    #D12052  magenta-red (positive extreme)

A yellow midpoint was tried first but the blue->yellow half runs through green,
which muddies the negative side and breaks the symmetry with the red half; a
white midpoint keeps a clean "no signal = blank" centre.

The :data:`thesis_cmap` / :data:`thesis_cmap_diverging` objects (and their
``_r`` reverses) and the :data:`THESIS_PALETTE` list are also importable
directly.  ``qcg`` / ``qcg_r`` remain registered as legacy aliases of
``thesis_cmap`` / ``thesis_cmap_r``.
"""

from __future__ import annotations

from matplotlib.colors import LinearSegmentedColormap

#: The palette colours, low value to high value.
THESIS_PALETTE = ["#03AED2", "#F8DE22", "#F45B26", "#D12052"]

# Sequential map.  The palette colours carry the bulk of the range; the ends are
# extended to a darker blue / deeper crimson so the absolute lowest and highest
# values stand out.
_SEQ_SEGMENTS = [
    (0.00, "#02566A"),  # darker shade of #03AED2 -- absolute low
    (0.13, "#03AED2"),
    (0.42, "#F8DE22"),
    (0.68, "#F45B26"),
    (0.89, "#D12052"),
    (1.00, "#9E0E3E"),  # deeper, more saturated crimson -- absolute high
]

# Diverging map: blue and red at the extremes, white as the neutral midpoint.
_DIV_SEGMENTS = [
    (0.0, "#03AED2"),
    (0.5, "#FFFFFF"),
    (1.0, "#D12052"),
]

#: Sequential colormap: dark blue (low) -> crimson (high).
thesis_cmap = LinearSegmentedColormap.from_list("thesis_cmap", _SEQ_SEGMENTS, N=256)
thesis_cmap_r = thesis_cmap.reversed()

#: Diverging colormap: blue (negative) -> white (zero) -> red (positive).
thesis_cmap_diverging = LinearSegmentedColormap.from_list(
    "thesis_cmap_diverging", _DIV_SEGMENTS, N=256
)
thesis_cmap_diverging_r = thesis_cmap_diverging.reversed()

# Legacy aliases (the sequential map used to be called "qcg").
_qcg = LinearSegmentedColormap.from_list("qcg", _SEQ_SEGMENTS, N=256)
_qcg_r = _qcg.reversed()

_ALL = (
    thesis_cmap,
    thesis_cmap_r,
    thesis_cmap_diverging,
    thesis_cmap_diverging_r,
    _qcg,
    _qcg_r,
)


def register(force: bool = True) -> None:
    """Register every thesis colormap with matplotlib's colormap registry.

    Called automatically on import.  Safe to call again; with ``force`` an
    already-registered name is replaced rather than raising.
    """
    import matplotlib

    registry = matplotlib.colormaps  # matplotlib >= 3.5
    for cm in _ALL:
        if cm.name in registry:
            if not force:
                continue
            registry.unregister(cm.name)
        registry.register(cm)


register()
