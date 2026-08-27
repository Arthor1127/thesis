#!/bin/bash
#$ -cwd
#$ -N amr_run1
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=128G
#$ -pe neworte 16
#$ -V

# Adaptive mesh refinement (AMR) test run - unlike the earlier probe
# scripts, this lets Palace actually SOLVE each AMR iteration (needed to
# compute the error indicator that drives refinement), so there is no
# kill-early trick here. Nobody will be watching this session live (the
# submitting machine is being turned off - see 2026-08-25 chat), so the
# whole thing is wrapped in a generous but FINITE `timeout` rather than
# relying on interactive supervision. Starts from an unbiased uniform
# 20um base mesh (~685k elements, confirmed cheap) and runs up to 3 AMR
# iterations with SaveAdaptMesh+SaveAdaptIterations=True, so the
# per-element `Indicator` field is available afterward in
# outputFiles/iteration{1,2,3}/ to see exactly where Palace refined.
#
# Memory: 128G, informed by evidence (the 2026-08-25 combined-6um-floor
# probe measured maxvmem=119GB for a much larger, uniformly-refined
# 38.6M-unknown problem; this AMR run starts far smaller (~685k-element
# base, Dorfler marking only refines a fraction of elements per
# iteration) so 128G is a safety margin, not an unfounded escalation.

set -uo pipefail

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
export PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"
export AMR_ITERS=3

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_test"
cd "$SCRIPT_DIR"

echo "=== Running on host: $(hostname) ==="
echo "=== NSLOTS: ${NSLOTS:-unset} ==="
echo "=== Start time: $(date) ==="

OUTDIR="$SCRIPT_DIR/amr_run1"
mkdir -p "$OUTDIR"

# 3-hour hard cap. If AMR hasn't finished by then, whatever iterations
# DID complete are still on disk (SaveAdaptIterations writes as it goes)
# for inspection later.
timeout 10800 python -u run_amr.py "$OUTDIR" 2>&1 | tee "$OUTDIR/amr.log"
STATUS=${PIPESTATUS[0]}

echo "=== End time: $(date) ==="
if [ "$STATUS" -eq 124 ]; then
    echo "!!! TIMED OUT after 3 hours - partial iteration data (if any) is still in $OUTDIR"
elif [ "$STATUS" -ne 0 ]; then
    echo "!!! FAILED with exit code $STATUS"
else
    echo "=== SUCCESS: AMR run completed within the time budget ==="
fi
