#!/bin/bash
#$ -cwd
#$ -N jjarr_fexport
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=48G
#$ -pe neworte 8
#$ -t 1-4
#$ -tc 4
#$ -V

# Sample each eigenmode's E-field onto a fixed 2D grid + chain-axis cut and save
# eig/field_data_n{N}.npz (read + sample with pyvista, no solve). One task per N.
#
# task:        1     2     3     4
#  n_junctions:10    23    40    100

set -uo pipefail
conda deactivate 2>/dev/null || true
conda activate qcg-quantum-design
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg PYVISTA_OFF_SCREEN=true

N_ARR=(10 23 40 100)
NJJ=${N_ARR[$((SGE_TASK_ID-1))]}

SCRIPT_DIR="$HOME/circuit_design/josephson_array_resonator"; cd "$SCRIPT_DIR"
echo "=== fexport task $SGE_TASK_ID: n=${NJJ}, host=$(hostname), $(date) ==="
python -u export_field_data.py \
    --run-dir "eig/n${NJJ}/eig_njj_${NJJ}_Lj_1.2nH" --z-um 2 2>&1
echo "=== task $SGE_TASK_ID done exit=$? $(date) ==="
