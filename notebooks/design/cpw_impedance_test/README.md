# CPW impedance/velocity validation

Validates `parametric_coupling.design_helper.DesignHelper.compute_cpw_geometric_parameters`
(the analytic CPW `Z0`, `v` formula) against SQDMetal/Palace FEM sims, as a
function of `cpw_width` and `cpw_gap` -- plus two supporting checks on the
method itself (length-linearity, film-thickness correction).

## Method

One straight CPW segment per grid point (`build_design.py`), but **not the
same geometry for both sims** -- `build_design.build_design(..., far_end=...)`
switches the far-end termination, and getting this backwards silently breaks
the capacitance extraction:

- **Capacitance** (`run_capacitance.py`, `far_end='open'`): `OpenToGround` at
  BOTH ends -- the trace is a single floating conductor, isolated from ground
  everywhere. Electrostatic trace<->ground Maxwell matrix -> total `C`,
  divided by the trace length -> `C'` (F/m). Using `far_end='short'` here
  would galvanically connect the trace to the ground plane (see
  `build_design.py`'s docstring and the geometry check below) -- trace and
  ground become ONE conductor, and there is no trace<->ground capacitance to
  extract; this is exactly what shipped in the first version of this pipeline
  and was caught before any cluster time was spent on it.
- **Self-inductance** (`run_inductance.py`, `far_end='short'`, the default):
  `OpenToGround` at the start, a REAL `ShortToGround` at the far end. A
  current is driven at the open end (via SQDMetal's U-clip, which locally
  shorts it to ground *for the sim only*) and returns through ground via the
  real short, so the closed loop encloses the CPW gap along the whole trace.
  Flux through that loop / current -> `L`, divided by the (margin-trimmed)
  length -> `L'` (H/m). This sim genuinely needs the galvanic short to close
  the loop -- it's the capacitance sim that must NOT have one.

Then `Z0 = sqrt(L'/C')`, `v = 1/sqrt(L'*C')`, compared point-by-point against
`build_design.analytic_reference(width, gap)`.

**Geometry sanity check** (re-run this after touching `build_design.py`):
```python
import shapely
from shapely.ops import unary_union
from build_design import build_design, _design_bounds_mm
d = build_design(20.0, 12.0, length_um=300.0, far_end='open')  # or 'short'
# buffer the "path" table (trace) and subtract everything else from the chip
# rect to get the actual ground polygon, then check trace.distance(ground):
# 0 for far_end='short' (must touch), == cpw_gap for far_end='open' (must not).
```

### Mesh resolution and AMR

Both sims size their base fine mesh at **half the smallest lateral feature**
(`0.5 * min(cpw_width, cpw_gap)`), then call `enable_mesh_refinement` (AMR,
on by default, see `AMR*` env vars / `--amr*` flags) to refine from there
toward a target relative error rather than trusting that fixed starting mesh
alone -- same `enable_mesh_refinement` API as
`squid_array_coupling/run_amr.py` (`num_iterations`, `max_DoFs`, `tolerance`,
`Dorfler_marking_fraction`). Both `run_capacitance.py` and `run_inductance.py`
additionally track `C`/`L` itself across AMR iterations (not just Palace's
internal FEM error estimator) and save that convergence trace into the result
JSON's `"amr"` key, mirroring `run_amr.py`'s per-iteration frequency check.

The open/short termination gaps are always included in the fine-mesh
component list -- omitting them collapses the isolating notch and silently
shorts an "open" end (the `sqdmetal-termination-gap-meshing` lesson).

### Zero-thickness metal vs. the film-thickness correction

The MAIN sweep (`gen_params.py` / `params.txt`) uses SQDMetal's default
zero-thickness ("quasi-3D") metal, matching
`scripts/python/island_cap/run_cap_island.py`'s finding that this is exact for
the in-plane field. `analytic_reference`'s finite-thickness correction terms
(Eqs. B6-B7 of `DesignHelper.compute_cpw_geometric_parameters`) are therefore
NOT exercised by the FEM side of the main sweep -- it fixes
`REF_THICKNESS_NM=100nm` on the analytic side only, as a representative value,
and treats the resulting small mismatch as part of the residual, not a bug.

The SEPARATE thickness sweep (`gen_thickness_params.py`,
`submit_*_thickness.sh`, `check_thickness_dependence.py`) instead runs
`FULL3D` (a real metal thickness) at thicknesses from 20 to 300nm, and
compares against `analytic_reference(..., thickness_nm=t)` at the SAME `t` --
an apples-to-apples test of that correction term specifically, at two
representative (width, gap) points rather than the whole grid (FULL3D meshing
is markedly more expensive).

### Finite-length fringe / length-linearity

Dividing one finite segment's total `C`/`L` by its length is only exact in the
limit length >> (cpw_width + 2*cpw_gap); at finite length each carries an
O(1) end-fringe correction that does not shrink with length in absolute terms,
only relative to it. The MAIN sweep's default `LENGTH_UM=3000` is chosen large
relative to the swept width/gap range to keep this a small systematic, but
that's an assumption, not a proof.

The SEPARATE length-linearity sweep (`gen_lengthcheck_params.py`,
`submit_*_lengthcheck.sh`, `check_length_linearity.py`) checks it directly:
runs the same (width, gap) point at several lengths and fits a straight line
through `C_total` (resp. `L_total`) vs. length. A good linear fit with a small
intercept (relative to a typical total) confirms `C_total/length` really is a
well-defined per-length quantity over the sweep's length range, rather than
something that happens to look consistent at one length. Uses two
representative points (this repo's usual 20/12um point, and the main grid's
`8/20um` corner -- the largest gap-to-feature ratio, where end-fringe
contamination should be largest) rather than the whole grid.

## Pipeline

1. **Local** (this machine, `qcg-quantum-design` or `quantum` conda env) --
   sanity-check geometry only, no FEM:
   ```
   python build_design.py            # smoke-builds one point, prints Z0/v
   python gen_params.py              # (re)writes params.txt if you change the grid
   python gen_lengthcheck_params.py  # writes params_lengthcheck.txt
   python gen_thickness_params.py    # writes params_thickness.txt
   ```
2. **Sync to the cluster** (never edit the cluster copy directly -- see
   `scripts/python/batch_test/sync_to_cluster.sh`'s own header for why):
   ```
   ../../../scripts/python/batch_test/sync_to_cluster.sh circuit_design/cpw_impedance_test --go
   ```
   Note: this only moves top-level `.py`/`.sh` files -- `params*.txt` are not
   `.py`/`.sh`, so either regenerate them on the cluster (`python gen_params.py`
   etc., inside a qsub job or, since these are pure-Python/no-FEM and
   near-instant, they're the one exception reasonable to run directly over SSH)
   or `rsync` them across explicitly.
3. **Cluster** -- all PALACE/FEM solves happen here, never on this machine
   (including "quick" validation runs -- a modest local eigenmode solve has
   previously driven this desktop into an unclean reboot), and package/env
   setup goes through a qsub job too, never run inline over the SSH frontend
   (see `submit_install_thesis_style.sh` for that pattern -- already run once;
   `thesis_style`/`thesis_cmaps` are now installed in the cluster's
   `qcg-quantum-design` env too, so plotting steps below can run there as an
   alternative to pulling results back). SSH in, `cd
   ~/circuit_design/cpw_impedance_test`, then, for each sweep, pilot a single
   task before launching the full array and resize `memoria_a_usar`/`h_rt`
   from what the pilot actually reports:
   ```
   qsub -t 1 submit_capacitance.sh              && qsub -t 1 submit_inductance.sh
   qsub -t 1 submit_capacitance_lengthcheck.sh   && qsub -t 1 submit_inductance_lengthcheck.sh
   qsub -t 1 submit_capacitance_thickness.sh     && qsub -t 1 submit_inductance_thickness.sh
   # then, once each pilot looks right:
   qsub submit_capacitance.sh              && qsub submit_inductance.sh
   qsub submit_capacitance_lengthcheck.sh  && qsub submit_inductance_lengthcheck.sh
   qsub submit_capacitance_thickness.sh    && qsub submit_inductance_thickness.sh
   ```
4. **Pull results back** (skip this if plotting on the cluster instead, per
   the note in step 3):
   ```
   for d in out_cap out_ind out_cap_lengthcheck out_ind_lengthcheck out_cap_thickness out_ind_thickness; do
     rsync -avz "arturo.ruiz@<cluster-host>:circuit_design/cpw_impedance_test/$d/results/" "$d/results/"
   done
   ```
5. **Analyze + plot**:
   ```
   python aggregate_and_plot.py           # main Z0(w,g)/v(w,g) validation
   python check_length_linearity.py       # C_total, L_total vs length -> per-length is well-defined
   python check_thickness_dependence.py   # Z0(t)/v(t) vs the analytic thickness correction
   ```
   All three use `thesis_style`/`thesis_cmaps` (`packages/qcg-quantum-design`,
   per `directives/figure-style.md`), which are pip-installed into this
   machine's own `quantum`/`qcg-quantum-design` conda envs AND (as of
   `submit_install_thesis_style.sh`) the cluster's `qcg-quantum-design` env --
   so these can run either here after pulling results back, or directly on
   the cluster. Pass `--no-plots` to just get the CSV (e.g. to sanity-check
   numbers right after a run without the plotting step).

## Files

| file | runs where | does |
|---|---|---|
| `build_design.py` | local + cluster | geometry (`far_end='open'|'short'`), integration area, analytic reference |
| `gen_params.py` / `gen_lengthcheck_params.py` / `gen_thickness_params.py` | local (or cluster, see step 2) | writes the three `params*.txt` grids |
| `run_capacitance.py` | cluster | one capacitance sim point (main / length-check / thickness sweep, via env vars) |
| `run_inductance.py` | cluster | one self-inductance sim point (ditto, via CLI flags) |
| `submit_capacitance*.sh` / `submit_inductance*.sh` | cluster (qsub) | SGE array jobs, one pair per sweep |
| `submit_install_thesis_style.sh` | cluster (qsub) | one-off: installs `thesis_style`/`thesis_cmaps` there |
| `aggregate_and_plot.py` | local or cluster | main sweep: merges results, validates vs. the analytic model, plots |
| `check_length_linearity.py` | local or cluster | length-check sweep: linearity + intercept fit, plots |
| `check_thickness_dependence.py` | local or cluster | thickness sweep: validates the B6-B7 correction term, plots |

Default main grid (`gen_params.py`): `cpw_width_um` in {8,12,16,20,26,34},
`cpw_gap_um` in {5,7,10,14,20} -- centred on this repo's common ~50 Ohm point
(`cpw_width=20um`, `cpw_gap=12um`, e.g. `notebooks/design/cap_coupler/`).
