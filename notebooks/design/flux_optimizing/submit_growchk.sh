#!/bin/bash
#$ -cwd
#$ -N flux_growchk
#$ -j y
#$ -S /bin/bash
#$ -q copahue,be_copahue,cpuINgpu
#$ -l memoria_a_usar=16G
#$ -pe neworte 8
#$ -t 1-36
#$ -tc 36
#$ -V

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

cd "$HOME/circuit_design/flux_optimizing"
OUT="growchk"
mkdir -p "$OUT"
K=$(printf '%02d' "$SGE_TASK_ID")
trap 'rm -rf "$HOME/circuit_design/flux_optimizing/$OUT/work_$K"' EXIT

echo "=== growchk $SGE_TASK_ID/36: host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 7200 python -u grow_check.py \
    --index "$SGE_TASK_ID" --outdir "$OUT" \
    --palace-bin "$PALACE_BIN" --num-cpus "$NSLOTS" \
    2>&1 | tee "$OUT/run_$K.log"
echo "=== task $SGE_TASK_ID exit=${PIPESTATUS[0]} $(date) ==="
