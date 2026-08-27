# Lessons Learned — Palace/SQDMetal on rocks7 (CentOS 7 / SGE)

Organized by symptom, so a future session can jump straight to the relevant
section instead of re-discovering these. See also: PALACE_INSTALL.md,
SQDMETAL_INSTALL.md, STATUS.md, setup_env.sh.

---

## Process discipline (read this first)

- **Claude has no direct access to the cluster.** It can only reason about
  what has been pasted into the conversation. Treat any claim about log
  content, memory figures, or results as unverified until the actual
  pasted output confirms it — this includes Claude's own prior messages
  in the same session. Once, mid-session, Claude stated specific mesh/
  memory numbers as observed fact without them having been pasted — they
  were fabricated predictions presented as results. Always paste real
  output; always ask for it if a claim seems to be getting ahead of the
  evidence.
- **A code fix Claude writes exists only in the chat until it is actually
  copied onto the cluster.** Multiple times a fix was written, discussed,
  and then a job was resubmitted using the *old* file because the copy
  step was skipped. Before resubmitting after any fix: `grep` the cluster
  file for the new logic (e.g. a new variable/arg name) to confirm the
  copy actually landed, not just that the job ran.
- **SGE job scripts run from a spooled copy, not the original file.**
  `$0` inside a `#$`-directive script resolves to a path under
  `/opt/gridengine/.../job_scripts/`, not wherever you actually saved the
  script. `cd "$(dirname "$0")"` silently breaks path-relative imports.
  Use `#$ -cwd` (runs the job from wherever `qsub` was invoked) and/or
  hardcode the real directory instead.

---

## Environment / conda / Spack

- **Conda must be fully deactivated before any Spack work.** Conda's own
  gcc/python/openssl/cmake can silently leak into Spack builds via `PATH`
  and contaminate the toolchain.
- **`spack load` and `conda activate` ordering matters.** `spack load`
  should come first, then `conda activate` — this way conda's `bin/`
  wins the `PATH` race for everyday tools (python, jupyter) while
  Spack's `mpirun` stays reachable further down `PATH`.
- **`unset LD_LIBRARY_PATH`** after both of the above, defensively, even
  though the real fix for the RPATH issue below was patching the
  binaries themselves — cheap insurance against it recurring.
- **Spack's default clone lands on `develop`, not a stable release.**
  Always `git checkout <latest vX.Y.Z tag>` explicitly.
- **Spack marks some patch versions deprecated** (hit this with
  `gcc@12.2.0`/`12.4.0`, and could recur with any package) — always
  `spack versions <pkg>` first and pick the highest non-deprecated one,
  or pass `--deprecated` to override.
- **Force a portable build target globally before building anything**:
  `spack config add packages:all:target:[x86_64]`. Without this, Spack
  auto-detects the *build node's* microarchitecture — a package built on
  a modern node can crash with "Illegal instruction" on this cluster's
  older Sandy-Bridge-era nodes.
- **CentOS 7 system externals (openssl, curl, etc.) often lack working
  pkg-config metadata.** Remove the broken `externals:` entries from
  `~/.spack/packages.yaml` and let Spack build its own from source
  instead of trusting the system version.
- **CMake's bootstrap can pick up a stray Qt5 config from an
  inherited/active conda environment via `PATH`**, even with
  `--no-qt-gui` passed, because the discovery happens through CMake's own
  test suite, not the main build. Fix: strip every `PATH` entry
  containing "conda" inside the job script before building CMake.

---

## RPATH poisoning (breaks Jupyter/zmq specifically)

**Symptom:** `import sqlite3` fails inside a Jupyter kernel with
`ImportError: .../libstdc++.so.6: version 'CXXABI_1.3.15' not found`,
while the identical `import sqlite3` works fine in a plain interactive
shell.

**Cause:** Several SQDMetal dependencies (scipy, pyproj, qutip, gdstk,
**pyzmq**) had to be built from source, using Spack's `gcc@12.5.0` as
`CC`/`CXX`. Their build systems baked an RPATH pointing at Spack's
`gcc-12.5.0/lib64` directly into the compiled `.so` files — RPATH takes
priority over `LD_LIBRARY_PATH` and cannot be overridden by any
environment variable. The first library with this bad RPATH to load in a
process (zmq loads very early in Jupyter kernel startup) "poisons" that
process's libstdc++ resolution for everything loaded afterward, even
unrelated libraries like `libicui18n` used by sqlite3.

**Fix:**
```bash
find $HOME/.conda/envs/qcg-quantum-design -iname "*.so*" -exec sh -c '
  readelf -d "$1" 2>/dev/null | grep -q "spack.*gcc-12.5.0" && patchelf --remove-rpath "$1" && echo "patched: $1"
' _ {} \;
```
Re-run periodically if new packages get compiled with the Spack
compiler in the future.

---

## Headless rendering (Qt / matplotlib)

- **Qt needs `QT_QPA_PLATFORM=offscreen`** on compute nodes (no display).
  Without it, `qiskit_metal`'s import can hard-abort the process
  (SIGABRT, not a catchable Python exception) — this crashed a Jupyter
  kernel silently with zero traceback until this was set.
- **`matplotlib.use("Agg", force=True)` must be called AFTER importing
  qiskit_metal/SQDMetal, not before.** qiskit_metal's own import chain
  force-switches the backend to `QtAgg` for its GUI; if Agg is set
  first, that later switch silently wins and overrides it.
- **VTK/PyVista 3D rendering does not work on these compute nodes at
  all** — no GPU, no X server, no EGL, no OSMesa. `pv.start_xvfb()` also
  doesn't exist in current PyVista. Do not spend time on this; use
  SQDMetal's own `PVDVTU_Viewer`/`retrieve_field_plots()` (matplotlib-
  based 2D slices) for inline field plots instead, or download ParaView
  output (`outputFiles/paraview/.../*.pvd`) to a local machine.
- **ParaView on a local machine may show N timesteps (one per MPI rank)
  instead of the correct number of cycles**, even though the same `.pvd`
  file reads correctly via PyVista/SQDMetal. Root cause not fully
  resolved — treat as a ParaView-side quirk, not a data problem, and
  fall back to `PVDVTU_Viewer` if it happens.

---

## Networking / proxy

- CNEA's proxy (`http://proxy.cnea.gob.ar:1280`) is needed both on the
  cluster (for pip/conda/git) and can interfere with SSH tunnel access
  to Jupyter from a local machine, since `localhost` traffic can get
  routed through the proxy if `no_proxy` isn't set — causing
  `ECONNREFUSED`/proxy-block errors even though the tunnel itself works
  (verify with `curl`, not just the browser/VSCode).
- VSCode's own proxy handling is separate from shell env vars — set
  `http.noProxy` and `http.proxySupport: "off"` explicitly in VSCode
  settings; shell-level `no_proxy` alone is not always sufficient for
  VSCode's Jupyter extension.
- Prefer the **direct compute-node URL**
  (`http://compute-X-Y.local:PORT/lab?token=...`) over an SSH tunnel when
  VSCode is already connected via Remote-SSH — it sidesteps the whole
  local-proxy problem since it never touches `localhost` on the laptop.

---

## SGE scheduling quirks

- `qstat -f -q <queue>` state flags (e.g. `E` for error) are the
  authoritative way to check if a specific node is actually healthy —
  more reliable than the "temporarily not available" text in
  `qstat -j <id>`'s scheduling info, which can be misleading/stale.
- A job that shows real progress in its log (mesh stats, etc.) and then
  reverts to `qw` almost certainly means the compute node it was on
  failed/was reclaimed and SGE requeued it — don't wait on it, kill and
  resubmit.
- **Priority aging**: a job's scheduling priority increases the longer
  it waits. Killing and resubmitting resets this. Only do it when there
  is real evidence the job is dead (see above) — not just because a
  single scheduling-info snapshot looks confusing while the target node
  is confirmed idle via `qhost`/`qstat -f`. A confirmed-idle node with a
  job stuck in `qw` for many hours is unusual and may be worth reporting
  to `soportefisica@cnea.gob.ar` rather than endlessly resubmitting.
- Parallel environment `neworte` has `allocation_rule=$pe_slots` — the
  **entire job's slots must come from a single host**. No multi-node MPI
  jobs are possible in this PE regardless of the 200-slot user quota;
  that quota is "up to 200 slots across many simultaneous jobs," not "one
  job can use 200 cores."
- `-l memoria_a_usar=N` is a **scheduling hint, not an enforced cap** —
  confirmed via `qconf -sq copahue` showing `s_vmem`/`h_vmem` = INFINITY.
  The actual kill on over-use comes from the Linux OOM-killer once real
  node RAM is exhausted, independent of what was requested.
- Real per-node RAM varies hugely in `copahue` (as of this session):
  `compute-4-7`/`4-8` ≈ 31G (small, avoid for anything beyond trivial
  jobs), most others 62–188G, `compute-4-18` ≈ 503G (largest). Check
  current values with `qhost | grep compute-4`, don't assume from memory
  — nodes can also enter an error state (`compute-4-9` had one).

---

## Mesh / geometry (SQDMetal + Palace)

- **`fine_mesh_components()` is whole-component-scoped, not
  feature-scoped.** Applying a fine `min_size` meant for a small feature
  (e.g. a 10um grounding strap) to the whole component it belongs to
  refines the *entire* component's edges — for a long meandered
  resonator, this caused a 20-50x element count blowup and repeated
  out-of-memory kills. Use `fine_mesh_in_rectangle(x1,y1,x2,y2,...)`
  (coordinates in **meters**, chip design uses **mm** internally) to
  target just the feature's real bounding box instead, extracted from
  `design.qgeometry.tables['poly']` after `design.rebuild()`.
- **A margin/padding around such a targeted box must never be
  clearance-derived from "distance to the same component's own trace
  geometry."** A grounding strap is physically attached to its trace by
  design — that distance is ~0 right at the attachment point, not a
  measure of clearance to a genuinely different neighboring feature.
  When this hit exactly 0.0, the resulting zero-width fine-mesh region
  was geometrically degenerate and crashed Palace almost instantly
  (signal 9, before any normal mesh-partitioning output even printed) —
  a mesh-validity crash, not an OOM, easy to misdiagnose as the same
  thing. **Use a small fixed margin instead** (3um worked well here).
- **A too-generous margin has the opposite, equally bad failure mode**:
  20um was enough to bridge into an actually-different neighboring
  meander strand (adjacent folds of the same trace can sit surprisingly
  close together), which made meshing *worse* (boundary elements ~1900 →
  ~42000) rather than better.
- **`fillet_resolution` (number of segments approximating curved
  bends) had ~zero measurable effect on element count** in this design
  (6 vs 12: 126,613,373 vs 126,606,114 elements) — not a useful memory
  lever, don't spend time tuning it.
- **The background mesh size (`min_size`/`max_size` on the whole domain)
  is the dominant memory lever**, far more impactful than the strap-
  specific refinement: going from 8um/100um to 15um/150um cut elements
  by 84% (126.6M → 19.7M) with no change to the strap's own resolution,
  and cut mesh-*generation* time (GMSH, before Palace/MPI even starts)
  by 6.3x (4362s → 687s).
- **Mesh generation and the actual FEM solve are separable and can be
  tested independently, fast**, since `gmsh` is a real importable Python
  module here (not a subprocess): call `driven_sim.prepare_simulation()`
  alone (no `.run()`), then query `gmsh.model.mesh.getNodes()` /
  `getElements()` directly for instant element/node counts — this avoids
  needing a full MPI/Palace launch (many minutes to hours) just to test
  a mesh-parameter change.
- **`solns_to_save`** controls how many frequency/mode steps get a full
  3D field snapshot written to disk for ParaView — it is a disk-I/O
  knob, not a memory lever for the solve itself. Don't set it to the
  full point count (e.g. 321) by default; a handful (5-10) is enough for
  visualization.
- **A driven S-parameter simulation is inherently much larger than the
  corresponding eigenmode simulation of just the resonator** — it
  additionally meshes the feedline and both launchpads, plus lumped-port
  boundary conditions, none of which the eigenmode-of-resonator-alone
  case touches. Don't expect them to have comparable memory footprints
  even with matching resonator-side settings.
- **`retrieve_data()`/`.run()` do not always propagate an MPI failure as
  a Python exception.** A `Killed` (signal 9, OOM) or a mesh-validity
  crash can leave `.run()` returning normally, with the actual failure
  only visible in the raw stdout/stderr log, and `retrieve_data()`
  returning `None` afterward. Both `build_and_eigenmode.py` and
  `build_and_driven.py` now explicitly check for this and fail loudly
  with a pointer to check for OOM/signal 9 in the log, rather than
  raising a confusing `NoneType` traceback.

---

## Solver order

- `solver_order=2` gives meaningfully better accuracy per element
  (important for eigenfrequency/Q precision, and for representing
  curved fillets) but costs roughly 4-8x more degrees of freedom than
  `solver_order=1` for the same mesh — a major, sometimes decisive,
  memory driver independent of mesh refinement choices.
