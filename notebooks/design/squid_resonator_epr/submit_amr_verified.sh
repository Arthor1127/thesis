#!/bin/bash
#$ -cwd
#$ -N amr_verified
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=64G
#$ -pe neworte 16
#$ -V

# AMR run built directly on top of the PROVEN local coarse_smoke_test
# mesh+config (Order=1, Target=7.0GHz, Tol=1e-5, MaxIts=300, 434131
# elements, converged in ~1 minute locally on 2026-08-24) - see
# 2026-08-25 chat: an earlier attempt to reconstruct this via a Python/
# SQDMetal script from a misremembered "20um uniform" mesh either
# segfaulted (Nonconformal=true) or failed to converge (Nonconformal=
# false), so this reuses the EXACT working mesh+json directly instead of
# regenerating it, with only a "Refinement" (AMR) block added on top.
# No Python/SQDMetal involved - straight Palace invocation.
#
# Nobody will be watching this session live (submitting machine is being
# turned off), so wrapped in a generous but finite timeout rather than
# relying on interactive supervision.

set -uo pipefail

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SIM_DIR="$HOME/circuit_design/parametric_coupling_test/amr_verified"
cd "$SIM_DIR"

echo "=== Running on host: $(hostname) ==="
echo "=== NSLOTS: ${NSLOTS:-unset} ==="
echo "=== Start time: $(date) ==="

timeout 10800 "$PALACE_BIN" -np "${NSLOTS:-16}" amr_verified.json 2>&1 | tee run.log
STATUS=${PIPESTATUS[0]}

echo "=== End time: $(date) ==="
if [ "$STATUS" -eq 124 ]; then
    echo "!!! TIMED OUT after 3 hours - partial AMR iteration data (if any) is still in outputFiles/iterationN/"
elif [ "$STATUS" -ne 0 ]; then
    echo "!!! FAILED with exit code $STATUS"
else
    echo "=== SUCCESS: AMR run completed within the time budget ==="
fi
