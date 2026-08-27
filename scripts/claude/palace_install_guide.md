# Installing Palace via Spack on rocks7 (CentOS 7 / SGE cluster)

This documents the full path to a working Palace install on the CAB/CNEA
`rocks7` cluster (frontend `10.73.25.223`, SGE scheduler, CentOS 7.9,
glibc 2.17, system gcc 4.8.5). Total elapsed effort involved several
false starts — this guide gives you the fixed sequence so you don't have
to rediscover each issue.

## 0. Before you start

- **Deactivate conda** for the whole process. Conda's own compilers,
  Python, and cmake will silently contaminate Spack builds via `PATH`,
  `CMAKE_PREFIX_PATH`, `LD_LIBRARY_PATH`, etc. `conda deactivate` alone
  is not always enough — see step 6.
- Confirm the CNEA proxy works from a compute node before submitting any
  long build (a stuck/misconfigured proxy wastes a full job cycle):
  ```bash
  cat > test_net.sh <<'EOF'
  #!/bin/bash
  #$ -cwd
  #$ -N test_net
  #$ -j y
  #$ -S /bin/bash
  #$ -q copahue
  #$ -l memoria_a_usar=1G
  #$ -V
  export http_proxy="http://proxy.cnea.gob.ar:1280"
  export https_proxy="http://proxy.cnea.gob.ar:1280"
  curl -sI https://github.com | head -5
  EOF
  qsub test_net.sh
  ```

## 1. Build a standalone Python for Spack itself

There is no system `python3` on this cluster (only `python2.7`), and
Spack needs Python ≥3.6 just to run. Build one from source with the
system gcc, isolated from conda:

```bash
mkdir -p ~/local/src ~/local/python3.9
cd ~/local/src
curl -O https://www.python.org/ftp/python/3.9.18/Python-3.9.18.tgz
tar xzf Python-3.9.18.tgz
cd Python-3.9.18
./configure --prefix=$HOME/local/python3.9
make -j$(nproc)
make install
```

Verify `ssl` and `zlib` work (Spack needs both for fetching):
```bash
~/local/python3.9/bin/python3 -c "import ssl, zlib; print('ok')"
```
(`sqlite3` may fail to import — irrelevant, ignore it.)

Point Spack at it explicitly and permanently:
```bash
echo 'export SPACK_PYTHON=$HOME/local/python3.9/bin/python3' >> ~/.bashrc
source ~/.bashrc
```

## 2. Install Spack itself, pinned to a stable release

```bash
git clone -c feature.manyFiles=true https://github.com/spack/spack.git ~/spack
cd ~/spack
git fetch --tags
git tag --list | sort -V | tail -10   # pick the latest stable vX.Y.Z tag
git checkout v1.2.2                   # use whatever is actually latest
echo 'source $HOME/spack/share/spack/setup-env.sh' >> ~/.bashrc
source ~/.bashrc
spack --version
```

Avoid `develop`/`main` — the default clone can land you there, which is
an unreleased, unstable branch. Always check out a tagged release.

## 3. Force a portable build target

The `copahue` queue spans wildly different hardware generations, from
2012 Sandy Bridge Xeons to 2020 Gold 6226R. If Spack auto-detects the
CPU of whichever node happens to build a package, that binary can crash
with "Illegal instruction" on older nodes. Fix this globally, once,
before building anything:

```bash
spack config add packages:all:target:[x86_64]
```

## 4. Bootstrap the concretizer

```bash
spack bootstrap now
```

May build `clingo` from source (CentOS 7 glibc is too old for prebuilt
bootstrap binaries) — this is normal.

## 5. Clean up broken external package registrations

Spack auto-detects system packages via `spack external find --all` and
registers them in `~/.spack/packages.yaml`. Several CentOS 7 system
packages are **too old / missing pkg-config metadata** and will break
downstream builds with errors like:
```
configure: error: openssl is a must but can not be found...
CMAKE_USE_SYSTEM_CURL is ON but a curl is not found!
```

The fix each time is the same: open `~/.spack/packages.yaml` and delete
the offending `externals:` block so Spack builds its own copy from
source instead of trusting the broken system one. We hit this for:
- `openssl` (two stale system versions, no working `.pc` file)
- `curl` (same issue, blocks CMake's `--system-curl` bootstrap)

If you hit a new "X is a must but cannot be found" / pkg-config error
during a build, check `~/.spack/packages.yaml` for an `externals:` entry
for that package pointing at `/usr`, and remove it.

## 6. Build a modern GCC, bootstrapped with the system compiler

```bash
cat > build_gcc.sh <<'EOF'
#!/bin/bash
#$ -cwd
#$ -N build_gcc
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=8G
#$ -pe mpi 8
#$ -V

unset PYTHONPATH PYTHONHOME
conda deactivate 2>/dev/null || true
export SPACK_PYTHON=$HOME/local/python3.9/bin/python3
export http_proxy="http://proxy.cnea.gob.ar:1280"
export https_proxy="http://proxy.cnea.gob.ar:1280"
export no_proxy="localhost,127.0.0.1,copahue,10.73.25.223"
source $HOME/spack/share/spack/setup-env.sh

spack install -j $NSLOTS gcc@12.5.0 %gcc@4.8.5
EOF
qsub build_gcc.sh
```

Notes:
- **Check `spack versions gcc` first.** Some patch versions (we hit this
  with 12.2.0 and 12.4.0) are marked deprecated in the current package
  repo and will fail immediately with "can only be satisfied by
  deprecated versions". Pick the highest non-deprecated version in the
  series you want (we landed on 12.5.0).
- Took ~30–35 minutes.

Register the new compiler once built:
```bash
spack compiler add $(spack location -i gcc@12.5.0 target=x86_64)/bin
spack compiler list
```

## 7. Build OpenMPI with the new compiler

```bash
cat > build_mpi.sh <<'EOF'
#!/bin/bash
#$ -cwd
#$ -N build_mpi
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=8G
#$ -pe mpi 8
#$ -V

unset PYTHONPATH PYTHONHOME
conda deactivate 2>/dev/null || true
export SPACK_PYTHON=$HOME/local/python3.9/bin/python3
export http_proxy="http://proxy.cnea.gob.ar:1280"
export https_proxy="http://proxy.cnea.gob.ar:1280"
export no_proxy="localhost,127.0.0.1,copahue,10.73.25.223"
source $HOME/spack/share/spack/setup-env.sh

spack install -j $NSLOTS openmpi@4.1.8 fabrics=none %gcc@12.5.0
EOF
qsub build_mpi.sh
```

Notes:
- **Spec ordering matters.** `%compiler` must come *after* variant flags
  like `fabrics=none`, not before, or Spack mis-parses which package the
  variant applies to.
- Check `spack versions openmpi` for deprecated versions, same as gcc.
- Took ~18 minutes (partly building its own OpenSSL after step 5's
  cleanup).

## 8. Build Palace itself

```bash
cat > build_palace.sh <<'EOF'
#!/bin/bash
#$ -cwd
#$ -N build_palace
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=16G
#$ -pe mpi 16
#$ -V

unset PYTHONPATH PYTHONHOME
conda deactivate 2>/dev/null || true
unset CMAKE_PREFIX_PATH PKG_CONFIG_PATH CPATH LIBRARY_PATH LD_LIBRARY_PATH

# Strip any conda-related entries from PATH — conda's own cmake/Qt5
# config files get picked up otherwise and break CMake's bootstrap
export PATH=$(echo "$PATH" | tr ':' '\n' | grep -v -i conda | paste -sd:)

export SPACK_PYTHON=$HOME/local/python3.9/bin/python3
export http_proxy="http://proxy.cnea.gob.ar:1280"
export https_proxy="http://proxy.cnea.gob.ar:1280"
export no_proxy="localhost,127.0.0.1,copahue,10.73.25.223"
source $HOME/spack/share/spack/setup-env.sh

spack install -j $NSLOTS palace@0.16.0 %gcc@12.5.0 ^openmpi@4.1.8
EOF
qsub build_palace.sh
```

Notes:
- Check `spack versions palace` and pick the latest non-deprecated
  release rather than `develop`.
- **The `unset CMAKE_PREFIX_PATH...` and conda PATH-stripping lines are
  essential.** Without them, CMake's own bootstrap process (triggered
  as a Palace dependency) picks up a stray Qt5 config file from an
  active/inherited conda environment (`.../miniconda3/lib/cmake/Qt5Gui`)
  and fails with a missing `GL/gl.h` error, even though CMake is invoked
  with `--no-qt-gui`. The Qt5 discovery happens through CMake's own test
  suite (`Tests/CMakeLists.txt`), not the main build, and isn't
  disabled by the GUI flag alone.
- Palace pulls in PETSc, SLEPc, MFEM, HYPRE, METIS, ParMETIS,
  SuperLU-dist, ScaLAPACK-equivalents, sundials, and CMake itself
  (bootstrapped from source). Once the compiler/MPI/openssl/curl layer
  is already built and cached, the full Palace chain took **~21
  minutes** end to end (a from-scratch run, e.g. on a fresh Spack
  instance, would be closer to 2-3 hours mostly dominated by the CMake
  bootstrap and PETSc build).

## 9. Verify and set up for regular use

```bash
spack find palace
spack load palace
palace --help
```

Generate environment modules so you don't need to source Spack's
setup-env.sh / spack load every session:
```bash
spack module tcl refresh
echo 'module use $HOME/spack/share/spack/modules/linux-centos7-x86_64' >> ~/.bashrc
source ~/.bashrc
module avail palace
module load palace
```

## Summary of root causes hit, for quick reference

| Symptom | Root cause | Fix |
|---|---|---|
| `spack: command not found` in fresh shell | setup-env.sh not sourced | add to `.bashrc` |
| gcc/openmpi "can only be satisfied by deprecated versions" | Spack repo marks some patch versions deprecated | `spack versions <pkg>`, pick highest non-deprecated |
| gcc binary "Illegal instruction" | Auto-detected build-node microarch (Haswell+) doesn't run on older cluster nodes (Sandy Bridge) | `spack config add packages:all:target:[x86_64]` before building anything |
| `libevent`/CMake "openssl/curl is a must but cannot be found" | Stale system openssl/curl registered as Spack externals, missing pkg-config `.pc` files | Remove the offending `externals:` block from `~/.spack/packages.yaml` |
| Spack "matches multiple packages" for gcc@12.5.0 | Old bootstrap-era gcc build with same version, different hash, still registered | Reference by hash: `spack location -i /hpbjuqh` |
| CMake bootstrap fails on missing `GL/gl.h` via Qt5Gui config, even with `--no-qt-gui` | Active/inherited conda environment's own cmake/Qt5 config leaking in via `PATH`, unrelated to `CMAKE_PREFIX_PATH` | Strip all conda entries from `PATH` inside the job script before sourcing Spack |
