#!/bin/bash
#$ -cwd
#$ -N flux_reconv
#$ -j y
#$ -S /bin/bash
#$ -q highmem,copahue,be_copahue
#$ -l memoria_a_usar=56G
#$ -pe neworte 12
#$ -t 1-25
#$ -tc 25
#$ -V

# Re-solve the 46 sweep points that hit SOLVER_MAXITS=1200 (all before the
# mid-sweep bump to 2000) without converging -- see reconverge_stragglers.py.
# One task = one straggler-bearing shard (25 shards, 1-4 points each), so
# concurrent tasks never append to the same CSV. Each point gets
# solver_maxits=6000, gmres_maxsize=300 (vs the sweep's 2000/250): only 46
# points total, so it's cheap to way overshoot the iteration budget. 56 G /
# 12 cores gives headroom over the ~40 G measured at MaxSize=250 -- MaxSize=300
# is ~20% more Krylov-vector memory on top of that.
#
# Submitted held on the two gap-fill arrays (-hold_jid) so it runs strictly
# last and its appended rows are the freshest -- aggregate_sweep.py's
# drop_duplicates(keep='last') then always prefers them regardless of which
# shard file they land in.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

cd "$HOME/circuit_design/flux_optimizing"
echo "=== flux_reconv task $SGE_TASK_ID/25: host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 36000 python -u reconverge_stragglers.py \
    --task-id "$SGE_TASK_ID" --n-tasks 25 \
    --outdir sweep \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    --solver-maxits 6000 --gmres-maxsize 300 \
    2>&1 | tee -a "sweep/reconv_$(printf '%02d' $SGE_TASK_ID).log"
STATUS=${PIPESTATUS[0]}
echo "=== flux_reconv task $SGE_TASK_ID exit=$STATUS $(date) ==="
exit $STATUS
