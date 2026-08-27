#!/bin/bash
#$ -cwd
#$ -N mesh_sweep
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=24G
#$ -pe neworte 2
#$ -t 1-56
#$ -tc 8
#$ -V

# One-at-a-time fine-mesh-size sweep, per region (readout / pad /
# squid_arm / junction), holding the other three regions at the coarse
# baseline (25um) while stepping the swept region down. Goal: find where
# GMSH mesh generation stops finishing within a 2-minute timeout, per
# region, so we know which region(s) actually drive the cost blowup
# (squid_arm was the suspect from the 2026-08-24 investigation) and how
# fine each can practically go.
#
# Modest resources per task (2 cores, 24G mem) since GMSH mesh generation
# is single-threaded and CPU/memory need scales with mesh size, which is
# capped hard by the 2-min timeout anyway - see the cluster-quota lesson
# from earlier in this session (don't request more than demonstrated need).
# -tc 8 caps concurrency at 8*2=16 slots, well under the 140 quota.

set -uo pipefail

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
export PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_test"
cd "$SCRIPT_DIR"

REGIONS=(readout pad squid_arm junction)
SIZES=(25 20 15 12 10 8 6 5 4 3 2.5 2 1.5 1)
N_SIZES=${#SIZES[@]}

IDX=$((SGE_TASK_ID - 1))
REGION_IDX=$((IDX / N_SIZES))
SIZE_IDX=$((IDX % N_SIZES))
REGION=${REGIONS[$REGION_IDX]}
SIZE=${SIZES[$SIZE_IDX]}

TASK_OUTDIR="$SCRIPT_DIR/sweep_run/${REGION}_${SIZE}um"
mkdir -p "$TASK_OUTDIR"

echo "=== Task $SGE_TASK_ID: region=$REGION min_size=${SIZE}um host=$(hostname) ==="

timeout 120 python -u mesh_sweep_step.py \
    --region "$REGION" --min-size-um "$SIZE" --outdir "$TASK_OUTDIR" \
    2>&1 | tee "$TASK_OUTDIR/task.log"
STATUS=${PIPESTATUS[0]}

if [ "$STATUS" -eq 124 ]; then
    echo "TIMEOUT after 120s"
    python3 -c "
import json
json.dump({'status': 'timeout', 'region': '$REGION', 'min_size_um': $SIZE}, open('$TASK_OUTDIR/result.json', 'w'), indent=2)
"
elif [ "$STATUS" -ne 0 ]; then
    echo "NONZERO_EXIT status=$STATUS"
fi

echo "=== Task $SGE_TASK_ID done, exit=$STATUS ==="
