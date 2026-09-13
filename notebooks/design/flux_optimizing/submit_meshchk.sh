#!/bin/bash
#$ -cwd
#$ -N flux_meshchk
#$ -j y
#$ -S /bin/bash
#$ -q copahue,be_copahue,cpuINgpu
#$ -l memoria_a_usar=16G
#$ -pe neworte 8
#$ -t 1-28
#$ -tc 28
#$ -V

# Mesh-convergence check for the flux sweep. 28 (point, mesh) tasks;
# row_<NN>.csv each. 8 cores/task so the fine (1.5 / 1.0 um) meshes finish
# in reasonable time.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

cd "$HOME/circuit_design/flux_optimizing"
OUT="meshchk"
mkdir -p "$OUT"

echo "=== meshchk task $SGE_TASK_ID/28: host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 10800 python -u mesh_check.py \
    --index "$SGE_TASK_ID" \
    --outdir "$OUT" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    2>&1 | tee "$OUT/run_$(printf '%02d' "$SGE_TASK_ID").log"
echo "=== task $SGE_TASK_ID exit=${PIPESTATUS[0]} $(date) ==="
