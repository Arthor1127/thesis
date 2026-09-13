#!/bin/bash
#$ -N mode_gap_corr
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=4G
#$ -pe neworte 1
#$ -o mode_gap_corr.log

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export MPLBACKEND=Agg
cd "$HOME/circuit_design/parametric_coupling_chip"
echo "=== mode_gap_corr host=$(hostname) $(date) ==="
python -u analysis/check_mode_gap_correlation.py 2>&1
ST=$?
echo "=== mode_gap_corr exit=$ST $(date) ==="
exit $ST
