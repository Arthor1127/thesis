#!/bin/bash
#$ -N plot_qc_fits
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=4G
#$ -pe neworte 1
#$ -o plot_qc_length_fits.log

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export MPLBACKEND=Agg
cd "$HOME/circuit_design/parametric_coupling_chip"
echo "=== plot_qc_fits host=$(hostname) $(date) ==="
python -u analysis/plot_qc_length_fits.py 2>&1
ST=$?
echo "=== plot_qc_fits exit=$ST $(date) ==="
exit $ST
