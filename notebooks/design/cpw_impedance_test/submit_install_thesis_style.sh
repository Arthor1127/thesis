#!/bin/bash
#$ -N thesis_style_install
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=2G
#$ -pe neworte 1
#$ -o install_thesis_style.log

# One-off: pip-installs the qcg-quantum-design package (thesis_style +
# thesis_cmaps, the house matplotlib style required by
# directives/figure-style.md) into the cluster's existing qcg-quantum-design
# conda env, so aggregate_and_plot.py-style figure generation can also run
# here instead of only on the local desktop. Source synced to
# ~/packages/qcg-quantum-design/ separately (rsync, not tracked by
# sync_to_cluster.sh, which only moves top-level .py/.sh).
#
# Per standing practice on this cluster, package installs are never run
# directly over SSH on the frontend node -- always via a qsub job, even for
# something this small.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh 2>/dev/null || true
conda activate qcg-quantum-design

echo "=== installing qcg-quantum-design (thesis_style/thesis_cmaps) into $(conda info --envs | grep '\*') ==="
pip install -e "$HOME/packages/qcg-quantum-design"
STATUS=$?

echo "=== verifying import ==="
python -c "import thesis_style, thesis_cmaps; print('OK:', thesis_style.__file__)"
VSTATUS=$?

echo "=== done: install_exit=$STATUS verify_exit=$VSTATUS $(date) ==="
exit $((STATUS + VSTATUS))
