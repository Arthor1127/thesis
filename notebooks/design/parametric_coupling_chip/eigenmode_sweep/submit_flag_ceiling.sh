#!/bin/bash
#$ -N flag_ceiling
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=4G
#$ -pe neworte 1
#$ -o flag_ceiling.log

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
conda activate qcg-quantum-design
cd "$HOME/circuit_design/parametric_coupling_chip"
echo "=== flag_ceiling host=$(hostname) $(date) ==="
python -u eigenmode_sweep/flag_eigenmode_ceiling.py 2>&1
ST=$?
echo "=== flag_ceiling exit=$ST $(date) ==="
exit $ST
