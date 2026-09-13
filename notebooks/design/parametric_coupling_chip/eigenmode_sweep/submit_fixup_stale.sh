#!/bin/bash
#$ -N qc_eig_fixup
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=4G
#$ -pe neworte 1
#$ -o fixup_stale.log

# JSON re-postprocessing only, no simulation -- still routed through qsub
# rather than the frontend, per standing convention (submit_refit.sh).

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
conda activate qcg-quantum-design
cd "$HOME/circuit_design/parametric_coupling_chip"
echo "=== qc_eig_fixup host=$(hostname) $(date) ==="
python -u eigenmode_sweep/fixup_stale_eigenmode_results.py 2>&1
ST=$?
echo "=== qc_eig_fixup exit=$ST $(date) ==="
exit $ST
