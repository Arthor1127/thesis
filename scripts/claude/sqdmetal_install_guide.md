# Installing SQDMetal on rocks7 (CentOS 7 / SGE cluster)

This documents the full path to a working SQDMetal install on the CAB/CNEA
`rocks7` cluster, on top of the Palace-via-Spack install (see the separate
Palace install guide). SQDMetal's own README only documents a Windows +
Anaconda install path — there is no Linux CI, no Docker, nothing tested
against an old-glibc HPC environment. This guide is the result of working
through every gap that leaves on CentOS 7.

**Prerequisite**: a working Palace install via Spack (`palace@0.16.0`, built
with `gcc@12.5.0`, `openmpi@4.1.8`) — see the Palace install guide. Also
assumes `conda`/miniconda3 is available on the cluster (`/share/apps/miniconda3`
here).

## 0. Clone the repo and create the conda environment

```bash
cd ~   # or wherever you want to keep it
git clone https://github.com/sqdlab/SQDMetal.git
conda create -n qcg-quantum-design python==3.11
conda activate qcg-quantum-design
```

Python 3.11 specifically — SQDMetal's `pyproject.toml` pins
`requires-python = ">=3.11,<3.12"`.

## 1. Remove the `mph` dependency (COMSOL support) before installing

SQDMetal supports two backends: COMSOL (via a Java bridge, `mph`/`jpype1`)
and Palace. If you don't have a COMSOL license on this cluster (you don't),
`mph`'s dependency chain will block the whole install for no benefit — see
step 5 below for why. Strip it out **before** running `pip install -e`:

```bash
cd ~/SQDMetal
grep -n "mph" pyproject.toml   # confirm the line number, e.g. "mph>=1.2.4",
cp pyproject.toml pyproject.toml.orig
sed -i '/"mph>=1.2.4",/d' pyproject.toml   # adjust the exact string if it differs
grep -n "mph" pyproject.toml   # should print nothing
```

## 2. Point the build at the Spack-built compiler

CentOS 7's system gcc (4.8.5) is too old to build several of SQDMetal's C
extensions (numpy, scipy, pyproj, etc.). Reuse the `gcc@12.5.0` you already
built for Palace:

```bash
GCC_BIN=$HOME/spack/opt/spack/linux-x86_64/gcc-12.5.0-hpbjuqhbk3gnanxoijhakmvhxfyxjv5g/bin
export CC=$GCC_BIN/gcc
export CXX=$GCC_BIN/g++
export PATH=$GCC_BIN:$PATH
$CC --version   # confirm: gcc 12.5.0, not the system 4.8.5
```

Adjust the hash in that path to whatever `spack find -l gcc@12.5.0` shows
on your system if it differs.

## 3. Pin numpy to a version with a real manylinux wheel

The default `numpy` PyPI resolves to a `linux_x86_64`-tagged wheel (no
`manylinux` compatibility tag), which builds/installs but is silently
broken for every downstream package that tries to build against it
(`numpy-config found: NO`, `Run-time dependency numpy found: NO`). This
breaks scipy, qutip, and anything else compiled against numpy's C API.

```bash
pip index versions numpy   # look for a version with a manylinux2014 wheel
pip install numpy==2.2.6   # confirmed manylinux_2_17/manylinux2014 wheel
```

Force this pin to propagate into pip's *isolated build* subprocesses too
(a plain `-c constraints.txt` flag does **not** reach nested build-dependency
installs — you need the environment variable):

```bash
echo "numpy==2.2.6" > /tmp/constraints.txt
export PIP_CONSTRAINT=/tmp/constraints.txt
```

## 4. Install system C-library dependencies via conda

Several of SQDMetal's dependencies need real C libraries that don't exist
on CentOS 7 and can't be `pip install`-ed as pure Python:

```bash
conda install -c conda-forge cmake meson meson-python ninja cython patchelf
conda install -c conda-forge proj gdal openblas vtk mesa-libgl-cos6-x86_64
conda install -c conda-forge openjdk=17
conda install -c conda-forge gmsh python-gmsh
```

Notes on each:
- **cmake**: needed by scipy, jpype1 (originally), and others; system has none.
- **proj/gdal**: needed by `pyproj` (a `geopandas` dependency, pulled in via
  `quantum-metal[full]`). `pyproj`'s build script needs `PROJ_DIR` set
  explicitly — see step 6.
- **openblas**: needed by scipy's meson build (`Dependency "OpenBLAS" not
  found`).
- **vtk**: needed by `pyvista`; PyPI's vtk wheel needs a newer glibc than
  CentOS 7 has, so let conda provide a compatible build.
- **openjdk=17**: `jpype1` (before its removal in step 1, and for whatever
  still needs it downstream) needs Java ≥9 to compile; system only has Java
  8 (`javac: invalid target release: 9`).
- **gmsh/python-gmsh**: the PyPI `gmsh` wheel is linked against a newer
  glibc than CentOS 7 provides (`GLIBC_2.23' not found`); conda's build
  targets the older glibc baseline.

## 5. Set PROJ_DIR

Conda may install this env under `~/.conda/envs/...` even if the base
conda binary lives under `/share/apps/miniconda3` — check which path is
actually correct for **your** env before setting this:

```bash
which proj   # confirms the real env path, e.g. ~/.conda/envs/qcg-quantum-design/bin/proj
export PROJ_DIR=$HOME/.conda/envs/qcg-quantum-design   # adjust to match `which proj`'s prefix
ls "$PROJ_DIR/include/proj.h"   # sanity check
```

Persist it:
```bash
echo 'export PROJ_DIR=$HOME/.conda/envs/qcg-quantum-design' >> ~/.bashrc
```

## 6. Install SQDMetal

```bash
cd ~/SQDMetal
pip install -e .
```

This pulls in the full dependency tree (pyvista, geopandas, scipy, qiskit-metal
via `quantum-metal[full]`, sphinx/nbsphinx doc tooling, etc.) — expect it to
take a while and print a lot of "Using cached ..." lines on a rerun since
most wheels/builds succeed and get cached after the first pass through the
fixes above.

If it still fails on something not covered here, the pattern that worked
throughout this whole install was: read the actual error (usually several
screens up from the final "ERROR: Failed to build X" line — meson/cmake
errors are verbose but tell you exactly what's missing), then either
`conda install -c conda-forge <thing>` if it's a system library, or `pip
install <thing>` if it's pure Python.

## 7. Fix remaining Python-level import failures

After `pip install -e .` succeeds, actually **importing** SQDMetal
(`import SQDMetal`) surfaces further gaps — these are runtime import
errors, not install failures, because `qiskit_metal` (a SQDMetal
dependency) eagerly imports every optional backend/GUI dependency at
module load time, regardless of which one you actually use.

**7a. Missing GUI/notebook packages** — `qiskit_metal` unconditionally
imports these; install as encountered:
```bash
pip install PySide6      # Qt backend qiskit_metal wants
pip install ipython
pip install pyEPR-quantum   # NOT "pyEPR" - that's a different, unrelated package
pip install pyaedt          # Ansys/AEDT interop - also unconditionally imported
```

**7b. jpype/mph imports scattered through SQDMetal's own COMSOL module**

Even after removing `mph` from `pyproject.toml` (step 1), `SQDMetal`'s own
source files under `SQDMetal/COMSOL/*.py` still directly `import mph` and
`import jpype...` at module level. Wrap every such import in try/except so
the (unused) COMSOL backend degrades gracefully instead of crashing the
whole `import SQDMetal`:

```bash
python3 - <<'EOF'
import re, glob

files = glob.glob("$HOME/SQDMetal/SQDMetal/COMSOL/*.py")
for p in files:
    s = open(p).read()
    orig = s
    s = re.sub(r'^import jpype\.types as jtypes$',
        'try:\n    import jpype.types as jtypes\nexcept ImportError:\n    jtypes = None',
        s, flags=re.MULTILINE)
    s = re.sub(r'^import jpype$',
        'try:\n    import jpype\nexcept ImportError:\n    jpype = None',
        s, flags=re.MULTILINE)
    s = re.sub(r'^import mph$',
        'try:\n    import mph\nexcept ImportError:\n    mph = None',
        s, flags=re.MULTILINE)
    if s != orig:
        open(p, "w").write(s)
        print("patched:", p)
EOF

find $HOME/SQDMetal/SQDMetal/COMSOL -name "__pycache__" -exec rm -rf {} +
```

**7c. Force matplotlib to a headless backend, and do it in the right order**

Compute nodes have no display. `qiskit_metal`'s own import chain
force-switches matplotlib to the interactive `QtAgg` backend (for its own
GUI) — if you call `matplotlib.use("Agg")` *before* importing
`qiskit_metal`, qiskit_metal's later internal switch silently overrides
you. The fix is to call it **after**:

```python
from SQDMetal.PALACE.Eigenmode_Simulation import PALACE_Eigenmode_Simulation
# ^ this import chain pulls in qiskit_metal, which force-switches to QtAgg

import matplotlib
matplotlib.use("Agg", force=True)   # must come AFTER the above import
```

Also set these at the job-script level as a backup, in case some other
code path (e.g. a subprocess Palace/GMSH spawns) tries to pick a GUI
backend independently:
```bash
export MPLBACKEND=Agg
export QT_QPA_PLATFORM=offscreen
```

## 8. Verify

```bash
python -c "import SQDMetal; print('SQDMetal imported OK')"
```

Should print cleanly with no traceback.

## 9. Load OpenMPI before actually running simulations

`import SQDMetal` succeeding does **not** mean simulations will run — Palace
itself needs `mpirun` on `PATH`, which isn't loaded by just activating the
conda env. Every session (interactive or in a job script) needs:

```bash
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
which mpirun   # confirm it resolves to the Spack-built OpenMPI, not empty
```

Without this, Palace fails fast and silently (`Error: Could not locate MPI
launcher`) rather than raising a clear Python exception — `driven_sim.run()`
returns in a fraction of a second having done nothing, and
`retrieve_data()` then fails on a missing output file.

## Summary of root causes hit, for quick reference

| Symptom | Root cause | Fix |
|---|---|---|
| `mph`/`jpype1` build chain fails (Java `-source 6` test target, ancient) | `mph` (COMSOL bridge) not needed for the Palace pathway | Remove from `pyproject.toml`; patch remaining `import mph`/`import jpype` sites in `SQDMetal/COMSOL/*.py` |
| `numpy-config found: NO` in scipy/qutip meson builds | Default numpy resolves to a non-manylinux-tagged wheel | Pin `numpy==2.2.6` (or any version with a real `manylinux2014` wheel) via `PIP_CONSTRAINT` env var, not just `-c` flag |
| `Dependency "OpenBLAS" not found` (scipy) | No BLAS/LAPACK on system | `conda install -c conda-forge openblas` |
| `ERROR: Invalid path for PROJ_DIR` / `proj.h` not found | conda env installed under `~/.conda/envs/`, not the base conda prefix you assumed | `export PROJ_DIR=$(dirname $(dirname $(which proj)))` (verify against `which proj`) |
| `javac: invalid target release: 9` | System Java is 1.8, jpype1 needs ≥9 | `conda install -c conda-forge openjdk=17`, set `JAVA_HOME` |
| `libm.so.6: version 'GLIBC_2.23' not found` (gmsh) | PyPI gmsh wheel needs newer glibc than CentOS 7's 2.17 | `conda install -c conda-forge gmsh python-gmsh` instead of pip |
| `ModuleNotFoundError: No module named 'ansys'/'IPython'/'PySide6'` on `import SQDMetal` | `qiskit_metal` eagerly imports every optional renderer/GUI backend | `pip install pyaedt ipython PySide6` |
| `ImportError: Cannot load backend 'QtAgg' ... 'headless' is currently running` | qiskit_metal force-switches matplotlib to Qt *after* your own `matplotlib.use("Agg")` call | Call `matplotlib.use("Agg", force=True)` **after** importing SQDMetal/qiskit_metal modules, not before |
| `Error: Could not locate MPI launcher` / driven-sim finishes in <1s with no output | `mpirun` not on `PATH` in the job/shell | `spack load openmpi@4.1.8 %gcc@12.5.0` before running any simulation |

## 10. Post-install fix required for Jupyter specifically (RPATH poisoning)

This gap was only discovered later, when running SQDMetal from a Jupyter
kernel rather than a plain script — worth doing proactively rather than
waiting to hit it.

**Symptom**: a Jupyter kernel crashes on `import sqlite3` (or similar)
with `libstdc++.so.6: version 'CXXABI_1.3.15' not found`, even though the
identical import works fine in a plain interactive shell.

**Cause**: `scipy`, `pyproj`, `qutip`, `gdstk`, and `pyzmq` were all built
from source (step 2 above pointed `CC`/`CXX` at Spack's `gcc@12.5.0`).
Their build systems baked an RPATH pointing at Spack's `gcc-12.5.0/lib64`
directly into the compiled `.so` files. RPATH takes priority over
`LD_LIBRARY_PATH` and cannot be overridden by any environment variable.
`pyzmq` loads very early in Jupyter kernel startup, so if it's the one
with the bad RPATH, it "poisons" libstdc++ resolution for the rest of
that process — even for unrelated libraries loaded later.

**Fix** — strip the bad RPATH from every affected file in the env:
```bash
conda install -c conda-forge patchelf
find $HOME/.conda/envs/qcg-quantum-design -iname "*.so*" -exec sh -c '
  readelf -d "$1" 2>/dev/null | grep -q "spack.*gcc-12.5.0" && patchelf --remove-rpath "$1" && echo "patched: $1"
' _ {} \;
```
Verify nothing is left:
```bash
find $HOME/.conda/envs/qcg-quantum-design -iname "*.so*" -exec sh -c 'readelf -d "$1" 2>/dev/null | grep -q "spack.*gcc-12.5.0" && echo "$1"' _ {} \;
```
(should print nothing). Re-run this check if any package gets rebuilt
from source with the Spack compiler in the future.

