#!/bin/bash
# Push this directory to the cluster and (optionally) submit a sweep array.
#   ./deploy.sh                              -- rsync only
#   ./deploy.sh submit <grid> <outdir> "q1,q2"
#     e.g. ./deploy.sh submit grid_omegar.json results_omegar "cpuINgpu,copahue,be_cpuINgpu,be_copahue"
#
# No --delete: this directory is shared with another agent's jobs.  Only the
# code/grid files are pushed; results_*/plots_*/logs stay put on both sides.
set -e
REMOTE=cnea-cluster:/home/arturo.ruiz/unitary/zak_phase
rsync -av \
  --exclude '_*' --exclude 'plots*' --exclude 'results*' \
  --exclude 'logs' --exclude 'sweep_grid.json' \
  --exclude '__pycache__' --exclude '*.pyc' \
  ./ "$REMOTE/"
echo "Synced -> $REMOTE"
if [ "$1" = "submit" ]; then
  ssh cnea-cluster "cd /home/arturo.ruiz/unitary/zak_phase && bash submit_zak_sweeps.sh '${2:-grid_omegar.json}' '${3:-results_omegar}' '${4:-cpuINgpu,copahue,be_cpuINgpu,be_copahue}'"
fi
