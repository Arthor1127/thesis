#!/bin/bash
# Push this directory to the cluster and optionally kick off the pipeline.
#   ./deploy.sh                 -- rsync only
#   ./deploy.sh setup           -- rsync + submit the conda-env build job
#   ./deploy.sh submit          -- rsync + submit the sweep array + plot job
#   ./deploy.sh all             -- rsync + setup, then submit (array holds on setup)
#
# No --delete: ~/unitary is shared with another agent's jobs.  Only code is
# pushed; results/ logs/ stay put on both sides.
set -e
REMOTE=cnea-cluster:/home/arturo.ruiz/unitary/cris
rsync -av \
  --exclude 'results' --exclude 'results_*' --exclude 'logs' --exclude 'plots*' \
  --exclude '__pycache__' --exclude '*.pyc' --exclude 'script-*' \
  ./ "$REMOTE/"
echo "Synced -> $REMOTE"

case "$1" in
  setup)  ssh cnea-cluster "cd /home/arturo.ruiz/unitary/cris && bash submit_setup.sh" ;;
  submit) ssh cnea-cluster "cd /home/arturo.ruiz/unitary/cris && bash submit_gaps.sh" ;;
  all)    ssh cnea-cluster "cd /home/arturo.ruiz/unitary/cris && bash submit_setup.sh && bash submit_gaps.sh" ;;
esac
