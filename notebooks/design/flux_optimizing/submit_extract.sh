#!/bin/bash
#$ -cwd
#$ -N flux_extract
#$ -j y
#$ -S /bin/bash
#$ -q highmem,copahue,be_copahue
#$ -l memoria_a_usar=8G
#$ -pe smp 2
#$ -V

# Post-process an already-finished eigenmode run's paraview output
# (extract_eigen_M.py: pyvista field read + integration, no Palace/FEM solve --
# still submitted as a job rather than run on the frontend, per house rule).
# TAG names the eigen/<TAG> directory (must already have result.json).

set -uo pipefail
conda deactivate 2>/dev/null || true
conda activate qcg-quantum-design
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg QISKIT_METAL_HEADLESS=1
PY=$HOME/.conda/envs/qcg-quantum-design/bin/python
cd "$HOME/circuit_design/flux_optimizing"
TAG="${TAG:?set TAG=<eigen subdir>}"
$PY -u extract_eigen_M.py "eigen/$TAG" ${CYCLE:+--cycle "$CYCLE"} 2>&1 | tee -a "eigen/$TAG/extract.log"
