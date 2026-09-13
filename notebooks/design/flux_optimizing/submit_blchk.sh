#!/bin/bash
#$ -cwd
#$ -N flux_blchk
#$ -j y
#$ -S /bin/bash
#$ -q copahue,be_copahue,cpuINgpu
#$ -l memoria_a_usar=48G
#$ -pe neworte 12
#$ -t 1-24
#$ -tc 24
#$ -V

# boundary-layer mesh convergence test. 4 points x 6 strategies. The ref05
# (uniform 0.5 um) tasks can hit ~10-20 M tets -> 12 cores + 48 G.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

cd "$HOME/circuit_design/flux_optimizing"
OUT="blchk"; mkdir -p "$OUT"
K=$(printf '%02d' "$SGE_TASK_ID")
trap 'rm -rf "$HOME/circuit_design/flux_optimizing/$OUT/work_$K"' EXIT

echo "=== blchk $SGE_TASK_ID/24: host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 14400 python -u bl_check.py \
    --index "$SGE_TASK_ID" --outdir "$OUT" \
    --palace-bin "$PALACE_BIN" --num-cpus "$NSLOTS" \
    2>&1 | tee "$OUT/run_$K.log"
echo "=== task $SGE_TASK_ID exit=${PIPESTATUS[0]} $(date) ==="
