#!/bin/bash
#$ -cwd
#$ -N check_mesh2
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=16G
#$ -pe neworte 2
#$ -V
set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
export PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"
cd "$HOME/circuit_design/parametric_coupling_test"
timeout 300 python -u check_base_mesh2.py runs_check2 2>&1
echo "=== exit: $? ==="
