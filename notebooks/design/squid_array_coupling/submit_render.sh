#!/bin/bash
#$ -cwd
#$ -N sqarr_render
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=96G
#$ -pe neworte 32
#$ -V

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
conda activate qcg-quantum-design
export MPLBACKEND=Agg QT_QPA_PLATFORM=offscreen

SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"
cd "$SCRIPT_DIR"
OF="$SCRIPT_DIR/amr_Lj600/amr_array/outputFiles"

python -u render_amr_refinement.py "$OF" "$SCRIPT_DIR/amr_Lj600/amr_refinement_volume.png" 2>&1
echo "=== render done, exit=$?, $(date) ==="
