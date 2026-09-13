#!/bin/bash
#$ -N qc_refit
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=4G
#$ -pe neworte 1
#$ -o refit.log

# Lightweight post-processing (curve fitting + plotting), NOT a FEM solve --
# still routed through qsub rather than run directly on the SSH frontend, per
# standing instruction (never run compute directly on the login node, even
# for small/quick things).

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export MPLBACKEND=Agg

cd "$HOME/circuit_design/parametric_coupling_chip"
echo "=== qc_refit host=$(hostname) $(date) ==="
python -u convergence_studies/refit_convergence.py 2>&1
ST=$?
echo "=== qc_refit exit=$ST $(date) ==="
exit $ST
