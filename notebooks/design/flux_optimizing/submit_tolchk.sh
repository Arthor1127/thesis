#!/bin/bash
#$ -cwd
#$ -N flux_tolchk
#$ -j y
#$ -S /bin/bash
#$ -q copahue,be_copahue,cpuINgpu
#$ -l memoria_a_usar=32G
#$ -pe neworte 12
#$ -t 1-8
#$ -tc 8
#$ -V
set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace"
cd "$HOME/circuit_design/flux_optimizing"; mkdir -p tolchk
K=$(printf '%02d' "$SGE_TASK_ID"); trap 'rm -rf "$HOME/circuit_design/flux_optimizing/tolchk/work_$K"' EXIT
echo "=== tolchk $SGE_TASK_ID: $(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 21600 python -u tol_check.py --index "$SGE_TASK_ID" --outdir tolchk --palace-bin "$PALACE_BIN" --num-cpus "$NSLOTS" 2>&1 | tee "tolchk/run_$K.log"
echo "=== exit=${PIPESTATUS[0]} $(date) ==="
