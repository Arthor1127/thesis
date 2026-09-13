#!/bin/bash
#$ -N qc_amr2
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=48G
#$ -pe neworte 16
#$ -o amr_pilot_noprom.log

# Isolation test for job 1938861's segfault (SIGSEGV, rank 0, ~92s into the
# solve, before the first AMR iteration finished): is it AMR+Driven crashing
# in general, or specifically AMR combined with OUR manually-patched PROM
# config (AdaptiveTol/AdaptiveMaxSamples, injected into the JSON after
# prepare_simulation() -- not part of SQDMetal's own tested API)? This run
# disables that patch (USE_PROM=0) and uses a small, narrow, DISCRETE sweep
# centered on the already-known fr=8.101GHz -- not trying to resolve a good
# Qc here, purely checking whether AMR+Driven survives at all without PROM.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_chip"
cd "$SCRIPT_DIR"
mkdir -p out_amr2

export COUPLING_SPACE_UM=14 COUPLING_LENGTH_UM=300
export FSTART_GHZ=8.05 FEND_GHZ=8.15 NPTS=41
export USE_PROM=0 SOLVER_TOL=1e-6 SOLVER_MAXITS=500
export AMR_ITERS=3 AMR_TOL=1e-2 AMR_MAX_DOFS=1500000 AMR_DORFLER=0.7
export TAG=amr_pilot_noprom PALACE_BIN OUTDIR="out_amr2"

echo "=== qc_amr2 (no PROM) host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u driven_sweep/run_qc_sim.py 2>&1
ST=$?
echo "=== qc_amr2 (no PROM) exit=$ST $(date) ==="
exit $ST
