#!/bin/bash
# SQDMetal / Palace environment setup.
#
# IMPORTANT: this must be SOURCED, not executed, since it needs to modify
# your current shell's environment (conda activate, exported variables).
# Running it as ./setup_env.sh in a subshell would have no effect on your
# actual terminal.
#
# Usage:
#   source ~/circuit_design/batch_test/setup_env.sh
# or, if you add the alias below to ~/.bashrc:
#   sqd_env

# --- Conda: source conda.sh first so `conda` is a shell function (activate/
# deactivate are no-ops otherwise in a non-interactive script). ---
source $HOME/miniconda3/etc/profile.d/conda.sh

# --- Start from a clean slate ---
conda deactivate 2>/dev/null || true
conda deactivate 2>/dev/null || true   # some setups nest more than one env

# --- Spack: gives us mpirun (OpenMPI 5.0.10, matching what Palace was built
# with - the actual spack root on this machine is ~/repo/spack, not ~/spack,
# and the installed openmpi is 5.0.10 %gcc@13.3.0, not 4.1.8 %gcc@12.5.0. ---
source $HOME/repo/spack/share/spack/setup-env.sh
spack load openmpi@5.0.10 %gcc@13.3.0

# --- Conda: must come AFTER spack load, so conda's python/bin/ wins the
# PATH race over Spack's paths for everyday tools, while mpirun (added by
# spack load above) stays reachable further down PATH. ---
conda activate qcg-quantum-design

# --- Unset LD_LIBRARY_PATH: several packages in this env (scipy, pyproj,
# qutip, gdstk, zmq) were originally built with Spack's gcc as CC/CXX and
# had its lib64 path baked into their RPATH, which could poison the whole
# process's libstdc++ resolution if LD_LIBRARY_PATH pointed there too.
# The RPATHs have since been stripped with patchelf, but this is cheap
# insurance against it recurring for anything rebuilt later. ---
unset LD_LIBRARY_PATH

# --- Headless rendering: no display on compute nodes. Qt needs the
# offscreen platform plugin or it can hard-abort (SIGABRT/SIGSEGV, not a
# catchable Python exception) as soon as qiskit_metal tries to initialize
# its GUI backend at import time. Matplotlib needs Agg for the same
# reason - and note qiskit_metal's own import silently force-switches the
# backend, so in scripts, matplotlib.use("Agg", force=True) must be called
# AFTER importing qiskit_metal/SQDMetal, not before. ---
export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg

# --- Palace binary path ---
export PALACE_BIN=$HOME/repo/spack/opt/spack/linux-zen3/palace-develop-blpbtxe43ne2or3yo7gou6ycmdwdxuzl/bin/palace

# --- Proxy (needed for pip/conda/anything reaching the internet) ---
export http_proxy="http://proxy.cnea.gob.ar:1280"
export https_proxy="http://proxy.cnea.gob.ar:1280"
export no_proxy="localhost,127.0.0.1,copahue,10.73.25.223"

# --- Sanity checks ---
echo "=================================================="
echo " SQDMetal / Palace environment check"
echo "=================================================="
echo "python:  $(which python 2>/dev/null)"
python --version 2>&1 | sed 's/^/         /'
echo "mpirun:  $(which mpirun 2>/dev/null)"
mpirun --version 2>&1 | head -1 | sed 's/^/         /'
echo "gcc:     $(which gcc 2>/dev/null)"
gcc --version 2>&1 | head -1 | sed 's/^/         /'
echo "PALACE_BIN: $PALACE_BIN"
if [ -x "$PALACE_BIN" ]; then
    echo "         (exists and is executable - OK)"
else
    echo "         !!! WARNING: not found or not executable !!!"
fi
echo "LD_LIBRARY_PATH: [${LD_LIBRARY_PATH:-<empty, good>}]"
echo "=================================================="

if [[ "$python" != *"qcg-quantum-design"* ]] && [[ "$(which python)" != *"qcg-quantum-design"* ]]; then
    echo "WARNING: python does not appear to be from the qcg-quantum-design env."
fi
