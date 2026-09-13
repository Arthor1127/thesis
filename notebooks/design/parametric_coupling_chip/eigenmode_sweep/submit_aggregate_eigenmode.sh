#!/bin/bash
#$ -N qc_eig_agg
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=4G
#$ -pe neworte 1
#$ -o aggregate_eigenmode.log

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export MPLBACKEND=Agg
cd "$HOME/circuit_design/parametric_coupling_chip"
echo "=== qc_eig_agg host=$(hostname) $(date) ==="
python -u eigenmode_sweep/aggregate_eigenmode_qc.py 2>&1
ST=$?
echo "=== qc_eig_agg exit=$ST $(date) ==="
exit $ST
