#!/bin/bash
#$ -cwd
#$ -N sqarr_smoke
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=48G
#$ -pe neworte 32
#$ -V

# First cluster run of the SQUID-array + readout-resonator EPR pipeline.
# Full n=23 geometry (so port-EPR.csv carries all 2*23 = 46 junction
# columns, the point of this design) but a deliberately COARSE mesh and
# only 4 modes - the goal is to confirm the pipeline runs end to end on
# the cluster and produces a well-formed epr_result.json with 46 inductive
# ports, not to get converged numbers. 48G is a guess; bump if it OOMs.
# Production mesh settings + memory will be set from what this run reports.

set -uo pipefail

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"
cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/smoke_n23"
mkdir -p "$OUTDIR"

echo "=== SQUID-array smoke: host=$(hostname), NSLOTS=$NSLOTS, $(date) ==="

timeout 3600 python -u build_and_eigenmode.py \
    --n-squids 23 \
    --L-j-nH 15 \
    --outdir "$OUTDIR" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    --number-of-freqs 4 \
    --starting-freq 3.0e9 \
    --bg-min-size-um 8 \
    --resonator-bg-min-size-um 8 \
    --array-bg-min-size-um 5 \
    --fillet-resolution 8 \
    2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}
echo "=== smoke done, exit=$STATUS, $(date) ==="
exit $STATUS
