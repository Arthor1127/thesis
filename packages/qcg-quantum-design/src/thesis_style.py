"""Thesis plotting style — palette, colour cycle, rcParams and helpers.

The colours are the same palette the :mod:`thesis_cmaps` colormaps are built
from.  The categorical **cycle** adapts it for line plots: the yellow is
deepened to a gold (pure ``#F8DE22`` is invisible as a thin line on white) and
two hues are added to reach six.

Quick start
-----------
::

    import thesis_style
    thesis_style.use()                     # standalone style, scienceplots fonts
    thesis_style.use(typography="newtx")   # Times-like maths via LaTeX
    thesis_style.use(typography=None)      # apply everything except fonts

    fig, axs = plt.subplots(2, 2, figsize=thesis_style.figsize(2, ratio=0.8))
    ...
    thesis_style.finalize(axs)
    thesis_style.panel_labels(axs)
    thesis_style.save(fig, "readout")      # -> readout.pdf + readout.png

What :func:`use` sets: the colour cycle, ``image.cmap = "thesis_cmap"``,
publication figure size / dpi, ticks pointing in with minor ticks and
top/right spines, a recessive grid, a frameless legend, mathtext tick
formatting, and ``savefig`` to a tight bbox.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

import thesis_cmaps  # noqa: F401  -- registers "thesis_cmap" etc.

# --------------------------------------------------------------------------
# Palette
# --------------------------------------------------------------------------

#: The six categorical colours, in cycle order (orange and red are kept apart
#: by the teal because they are only ~ΔE 14 for normal colour vision).
CYCLE = ["#03AED2", "#F45B26", "#12A594", "#D12052", "#E0A500", "#6A3D9A"]
CYCLE_NAMES = ["blue", "orange", "teal", "red", "gold", "purple"]

#: Linestyles and markers paired with :data:`CYCLE`, so series stay distinct in
#: greyscale / for colour-vision deficiency (used by ``use(bw=True)``).
CYCLE_LINESTYLES = ["-", "--", "-.", ":", "-", "--"]
CYCLE_MARKERS = ["o", "s", "^", "D", "v", "P"]

#: Named access, e.g. ``thesis_style.COLORS["teal"]``.
COLORS = dict(zip(CYCLE_NAMES, CYCLE))

#: ``ListedColormap`` of the cycle, for categorical ``imshow`` / bar charts.
thesis_cat_cmap = ListedColormap(CYCLE, name="thesis_cat")


def color_cycler(bw: bool = False):
    """Return the :class:`~cycler.Cycler` used for ``axes.prop_cycle``.

    With ``bw`` the cycle also varies linestyle and marker so the series are
    distinguishable without colour.
    """
    cyc = mpl.cycler(color=CYCLE)
    if bw:
        cyc = cyc + mpl.cycler(linestyle=CYCLE_LINESTYLES) + mpl.cycler(marker=CYCLE_MARKERS)
    return cyc


# --------------------------------------------------------------------------
# Figure sizes  (APS single / 1.5 / double column, inches)
# --------------------------------------------------------------------------

_COLUMN_WIDTH = {1: 3.375, 1.5: 5.0, 2: 6.75}


def figsize(cols: float = 1, ratio: float = 0.68, height: float | None = None):
    """``(width, height)`` in inches for a ``cols``-column figure.

    ``ratio`` is height/width; pass ``height`` to override it directly.
    """
    w = _COLUMN_WIDTH.get(cols, cols)
    return (w, height if height is not None else w * ratio)


# --------------------------------------------------------------------------
# Typography
# --------------------------------------------------------------------------

_USETEX_PREAMBLE = {
    "newtx": r"\usepackage{amsmath}\usepackage{newtxtext,newtxmath}",
    "txfonts": r"\usepackage{amsmath}\usepackage{txfonts}",
    "lmodern": r"\usepackage{amsmath}\usepackage{lmodern}",
    "palatino": r"\usepackage{amsmath}\usepackage{mathpazo}",
    "cm": r"\usepackage{amsmath}",
}
_MATHTEXT = {
    "stix": ("serif", "stix", ["STIX Two Text", "STIXGeneral", "DejaVu Serif"]),
    "dejavu-serif": ("serif", "dejavuserif", ["DejaVu Serif"]),
    "sans": ("sans-serif", "stixsans", None),
    "dejavu-sans": ("sans-serif", "dejavusans", None),
    # exactly what `plt.style.use(["science", "no-latex"])` sets: STIXGeneral
    # text with Computer-Modern-style mathtext, no LaTeX.
    "scienceplots": ("serif", "cm", ["STIXGeneral", "DejaVu Serif"]),
}

#: Typography keys accepted by :func:`use`.
TYPOGRAPHIES = sorted(set(_USETEX_PREAMBLE) | set(_MATHTEXT))


def _font_rc(typography: str | None) -> dict:
    if typography is None:
        return {}
    if typography in _USETEX_PREAMBLE:
        fam = "sans-serif" if typography == "sans" else "serif"
        return {
            "text.usetex": True,
            "font.family": fam,
            "text.latex.preamble": _USETEX_PREAMBLE[typography],
        }
    if typography in _MATHTEXT:
        family, fontset, serif = _MATHTEXT[typography]
        rc = {"text.usetex": False, "font.family": family,
              "mathtext.fontset": fontset, "axes.formatter.use_mathtext": True}
        if serif:
            rc["font.serif"] = serif
        return rc
    raise ValueError(f"unknown typography {typography!r}; choose from {TYPOGRAPHIES} or None")


# --------------------------------------------------------------------------
# The style
# --------------------------------------------------------------------------

def rc(typography: str | None = "scienceplots", *, bw: bool = False, grid: bool = True,
       cols: float = 1, ratio: float = 0.68) -> dict:
    """Build (but do not apply) the rcParams dict for the thesis style."""
    d = {
        "figure.figsize": figsize(cols, ratio),
        "figure.dpi": 120,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.transparent": False,
        "savefig.format": "pdf",

        "axes.prop_cycle": color_cycler(bw),
        "image.cmap": "thesis_cmap",
        "axes.linewidth": 0.8,
        "axes.grid": grid,
        "grid.color": "0.85",
        "grid.linewidth": 0.6,
        "grid.alpha": 0.7,
        "axes.axisbelow": True,

        "lines.linewidth": 1.8,
        "lines.markersize": 5.0,
        "lines.markeredgewidth": 0.8,
        "errorbar.capsize": 2.5,
        "scatter.edgecolors": "none",

        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,

        "legend.frameon": False,
        "legend.handlelength": 1.6,
        "legend.columnspacing": 1.2,

        "font.size": 9,
        "axes.titlesize": 9,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
    }
    d.update(_font_rc(typography))
    return d


def use(typography: str | None = "scienceplots", *, bw: bool = False, grid: bool = True,
        cols: float = 1, ratio: float = 0.68) -> None:
    """Apply the thesis style to the global rcParams.

    Parameters
    ----------
    typography : one of :data:`TYPOGRAPHIES` or ``None``.
        ``"scienceplots"`` (default) reproduces ``["science", "no-latex"]``:
        STIXGeneral text with Computer-Modern mathtext, no LaTeX.  ``"stix"``,
        ``"dejavu-serif"``, ``"sans"`` and ``"dejavu-sans"`` are the other
        no-LaTeX options.  ``"newtx"`` / ``"txfonts"`` / ``"lmodern"`` /
        ``"palatino"`` / ``"cm"`` render maths through LaTeX (needs a TeX
        install); ``"cm"`` is what plain ``["science"]`` uses.  ``None`` leaves
        the fonts untouched.
    bw : also cycle linestyle + marker (greyscale / CVD safe).
    grid : draw the recessive grid.
    """
    if "thesis_cat" not in plt.colormaps():
        mpl.colormaps.register(thesis_cat_cmap)
    plt.rcParams.update(rc(typography, bw=bw, grid=grid, cols=cols, ratio=ratio))


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

def _flat(axs):
    if isinstance(axs, plt.Axes):
        return [axs]
    return list(getattr(axs, "flat", axs))


def finalize(axs) -> None:
    """Tidy one or many axes: minor ticks on, ticks on all four sides."""
    for ax in _flat(axs):
        ax.minorticks_on()
        ax.tick_params(which="both", top=True, right=True)


def panel_labels(axs, labels=None, *, loc="upper left", pad=0.04, weight="bold",
                 fontsize=None, **text_kw):
    """Stamp ``(a) (b) (c) …`` on a set of axes (in axes-fraction coords)."""
    axl = _flat(axs)
    if labels is None:
        labels = [f"({chr(ord('a') + i)})" for i in range(len(axl))]
    x, ha = (pad, "left") if "left" in loc else (1 - pad, "right")
    y, va = (1 - pad, "top") if "upper" in loc else (pad, "bottom")
    for ax, lab in zip(axl, labels):
        ax.text(x, y, lab, transform=ax.transAxes, ha=ha, va=va,
                weight=weight, fontsize=fontsize, **text_kw)


def save(fig, name, *, formats=("pdf", "png"), dpi=600, directory="."):
    """Write ``fig`` to ``name`` once per extension in ``formats``."""
    import os

    paths = []
    for ext in formats:
        p = os.path.join(directory, f"{name}.{ext}")
        fig.savefig(p, dpi=dpi)
        paths.append(p)
    return paths


def highlight(ax, index, *, dim="0.78", dim_lw=1.0):
    """Grey out every line on ``ax`` except ``index``, which is brought to front."""
    lines = ax.get_lines()
    for i, ln in enumerate(lines):
        if i == index:
            ln.set_zorder(max(l.get_zorder() for l in lines) + 1)
        else:
            ln.set_color(dim)
            ln.set_linewidth(dim_lw)
            ln.set_zorder(ln.get_zorder() - 1)
