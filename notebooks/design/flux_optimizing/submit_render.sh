#!/bin/bash
#$ -cwd
#$ -N flux_render
#$ -j y
#$ -S /bin/bash
#$ -q cpuINgpu,be_cpuINgpu,copahue,be_copahue
#$ -l memoria_a_usar=3G
#$ -pe smp 2
#$ -t 1-100
#$ -tc 80
#$ -V
# Geometry-only render (no Palace). 1-slot-ish array tasks fit in any gap.
set -uo pipefail
conda deactivate 2>/dev/null || true
conda activate qcg-quantum-design
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg QISKIT_METAL_HEADLESS=1
PY=$HOME/.conda/envs/qcg-quantum-design/bin/python
cd "$HOME/circuit_design/flux_optimizing"
$PY -u render_geometry.py --mode tiles --outdir geometry_gallery --nprocs "$NSLOTS" --shard "$SGE_TASK_ID" --nshards 100 2>&1
