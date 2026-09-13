# qcg-quantum-design

Single source of truth for the group's `PHOTONICChain` (SSH photonic chain of
transmons with capacitive + cavity-mediated inductive coupling and a cavity)
and its second-quantization / Bogoliubov engine.

Top-level modules exposed:

- `ssh_chain_sq` — `PHOTONICChain`, `SecondQuantizationEngine`, `CompositeEngine`,
  `TensorTerm`, `Superposition`, `QuantumOperator`, `compute_position_exponential_sparse`.
- `sq_engine` — compatibility shim re-exporting the engine half of `ssh_chain_sq`.
- `thesis_cmaps` — the `thesis_cmap` / `thesis_cmap_diverging` matplotlib
  colormaps (see below). `qcg_cmaps` stays importable as a legacy shim.
- `thesis_style` — categorical colour cycle, standalone rcParams style
  (`thesis_style.use()`) and plotting helpers (see below).

Full API reference (every class and method, with physical scope):
[`docs/reference.md`](docs/reference.md).

## Install (editable, into whichever conda env)

```bash
conda run -n <env> pip install -e /home/ruiz/Documents/thesis/packages/qcg-quantum-design
```

Currently installed editable into `quantum` and `qcg-quantum-design`.

CuPy is used automatically when importable (`pip install -e '.[gpu]'`);
otherwise everything runs on NumPy.

## Colormaps

```python
import thesis_cmaps          # registers the colormaps with matplotlib
plt.imshow(data, cmap="thesis_cmap")            # sequential
plt.imshow(signed, cmap="thesis_cmap_diverging",
           norm=TwoSlopeNorm(vcenter=0))        # diverging
```

Both are built from the group palette
`#03AED2` blue · `#F8DE22` yellow · `#F45B26` orange · `#D12052` red.

- **`thesis_cmap`** — sequential, low → high. The palette runs blue → yellow →
  orange → red, with the ends extended to a darker blue (`#02566A`) and a deeper
  crimson (`#9E0E3E`) so the extremes stand out.
- **`thesis_cmap_diverging`** — diverging, `#03AED2` blue (−) → white (0) →
  `#D12052` red (+). A yellow midpoint was tried first but the blue→yellow half
  runs through green; white keeps the two halves symmetric.

Each has a `_r` reverse. `thesis_cmaps.thesis_cmap`,
`thesis_cmaps.thesis_cmap_diverging` (and `_r`), `thesis_cmaps.THESIS_PALETTE`
and `thesis_cmaps.register()` are importable directly. `qcg` / `qcg_r` stay
registered as legacy aliases of `thesis_cmap`; `import qcg_cmaps` still works.

## Plot style

```python
import thesis_style as ts
ts.use()                        # standalone rcParams style; scienceplots fonts
ts.use(typography="newtx")      # Times-like maths via LaTeX (needs a TeX install)

fig, axs = plt.subplots(2, 2, figsize=ts.figsize(2, ratio=0.8))
...
ts.finalize(axs)                # minor ticks on, ticks all four sides
ts.panel_labels(axs)            # (a) (b) (c) ...
ts.save(fig, "readout")         # -> readout.pdf + readout.png
```

- **Colour cycle** — `#03AED2` blue · `#F45B26` orange · `#12A594` teal ·
  `#D12052` red · `#E0A500` gold · `#6A3D9A` purple. The colormap palette adapted
  for lines: the yellow is deepened to a line-safe gold, two hues added, orange
  and red kept apart. `ts.CYCLE` / `ts.COLORS` / `ts.CYCLE_NAMES`; `thesis_cat` is
  the matching `ListedColormap`. `ts.use(bw=True)` also cycles linestyle + marker.
- **`ts.use()`** sets a *standalone* style (not layered on scienceplots): the
  cycle, `image.cmap="thesis_cmap"`, APS column figure sizes, 600-dpi pdf saves
  with a tight bbox, ticks-in with minor ticks and top/right spines, a recessive
  grid, a frameless legend, and mathtext tick labels. `typography=` picks the
  font: `scienceplots` (default, == `["science", "no-latex"]`: STIXGeneral text
  + Computer-Modern mathtext, no LaTeX), `stix`, `dejavu-serif`, `sans`, or the
  LaTeX options `newtx` / `txfonts` / `lmodern` / `palatino` / `cm` (`cm` ==
  plain `["science"]`); `None` leaves fonts alone.
- **Helpers** — `figsize`, `finalize`, `panel_labels`, `highlight(ax, i)`, `save`.

A showcase notebook covering the colormaps and the style lives at
[`notebooks/unitary/thesis_colormap_showcase.ipynb`](../../notebooks/unitary/thesis_colormap_showcase.ipynb).

## Note on local shadowing

`notebooks/exploration/unitary/` keeps its own older `ssh_chain_sq.py` /
`sq_engine.py`; code run with that directory as the working dir imports those
local copies, not this package.
