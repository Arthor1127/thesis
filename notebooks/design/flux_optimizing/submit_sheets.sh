#!/bin/bash
#$ -cwd
#$ -N flux_sheets
#$ -j y
#$ -S /bin/bash
#$ -q cpuINgpu,be_cpuINgpu,copahue,be_copahue
#$ -l memoria_a_usar=6G
#$ -pe smp 1
#$ -V
# Compose the 16 contact sheets from the rendered tiles + merge SUSPECT.txt.
# Submit held on the tile array:  qsub -hold_jid flux_render submit_sheets.sh
set -uo pipefail
conda deactivate 2>/dev/null || true
conda activate qcg-quantum-design
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg QISKIT_METAL_HEADLESS=1
PY=$HOME/.conda/envs/qcg-quantum-design/bin/python
cd "$HOME/circuit_design/flux_optimizing"
$PY -u render_geometry.py --mode sheets --outdir geometry_gallery --tile-w 1400 2>&1
echo "tiles rendered: $(find geometry_gallery/tiles -name '*.png' | wc -l) / 2400"
touch geometry_gallery/SHEETS_DONE.txt
