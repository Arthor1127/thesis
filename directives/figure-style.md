# Directive — make every figure with `thesis_style`

**All figures for this thesis are produced with the `thesis_style` /
`thesis_cmaps` modules from the `qcg-quantum-design` package.** Do not hand-pick
colours, do not use `jet` / `rainbow` / `viridis`-by-default, and do not layer
`scienceplots` on top. Future sessions: apply this to any plotting task for this
project unless the user explicitly asks otherwise.

## Where it lives

`packages/qcg-quantum-design/src/thesis_style.py` and `thesis_cmaps.py`. The
package is `pip install -e`'d into **both** the `quantum` and
`qcg-quantum-design` conda envs and points at the repo copy, so editing the repo
file updates both envs. A running Jupyter kernel must be **restarted** after a
fresh `pip install -e` or after a new module is added to the package.

Full walk-through: `notebooks/unitary/thesis_colormap_showcase.ipynb`.
API summary: `packages/qcg-quantum-design/README.md` and `docs/reference.md`.

## How to use it

```python
import numpy as np
import matplotlib.pyplot as plt
import thesis_style as ts

ts.use()                         # apply the standalone style (call once, early)

fig, axs = plt.subplots(2, 2, figsize=ts.figsize(2, ratio=0.8))
# ... plot ...
ts.finalize(axs)                 # minor ticks on, ticks on all four sides
ts.panel_labels(axs)             # (a) (b) (c) ...
ts.save(fig, "readout")          # writes readout.pdf + readout.png
```

`ts.use()` is a **standalone** rcParams style (not built on `scienceplots`). It
sets: the categorical colour cycle, `image.cmap = "thesis_cmap"`, APS
single-column figure size, 600-dpi PDF saves with a tight bbox, ticks pointing
in with minor ticks and top/right spines, a light grid behind the data, a
frameless legend, and mathtext tick labels.

### Colours

- **Line / categorical data** — just use the default cycle (blue, orange, teal,
  red, gold, purple). Named access: `ts.COLORS["teal"]`, list: `ts.CYCLE`.
  For greyscale / colour-blind safety call `ts.use(bw=True)` (also cycles
  linestyle + marker).
- **Sequential fields** (`imshow`, `pcolormesh`, `contourf`) — `cmap="thesis_cmap"`
  (the default), low → high.
- **Signed fields** where zero matters — `cmap="thesis_cmap_diverging"` with a
  symmetric norm, e.g. `norm=TwoSlopeNorm(vcenter=0, vmin=-a, vmax=a)`.
- **Categorical images / grouped bars** — `cmap="thesis_cat"` or `color=ts.CYCLE`.

The palette is `#03AED2` blue · `#F8DE22` yellow · `#F45B26` orange · `#D12052`
red; the cycle deepens the yellow to a line-safe gold `#E0A500` and adds teal +
purple. `thesis_cmaps.THESIS_PALETTE` holds the four originals.

### Typography

Default `typography="scienceplots"` reproduces
`plt.style.use(["science", "no-latex"])` — STIXGeneral text with Computer-Modern
mathtext, no LaTeX. Override per call if a specific target needs it:
`ts.use(typography="newtx")` (Times-like maths, needs LaTeX), or `"cm"`,
`"stix"`, `"sans"`, `"lmodern"`, `"palatino"`, `"txfonts"`, `"dejavu-serif"`,
`"dejavu-sans"`, or `None` to leave fonts untouched.

### Helpers

| call | does |
|------|------|
| `ts.figsize(cols, ratio=0.68)` | `(w, h)` inches; `cols` = 1 / 1.5 / 2 APS columns |
| `ts.finalize(axs)` | minor ticks on, ticks on all four sides (one or many axes) |
| `ts.panel_labels(axs)` | stamp `(a) (b) (c) …` |
| `ts.highlight(ax, i)` | grey every line but the `i`-th, bring it to the front |
| `ts.save(fig, name)` | write `name.pdf` and `name.png` in one call |

## If `thesis_style` is missing from an env

```bash
conda run -n <env> pip install -e /home/ruiz/Documents/thesis/packages/qcg-quantum-design
```

then restart the kernel.
