#!/bin/bash
# Submits the entire sweep pipeline in one shot: Stage A (eigenmode,
# 10 tasks), Stage B (driven S21, 10 tasks - one continuous 4.5-20.5 GHz
# sweep per position), then Stage C (aggregation + final plots), which
# automatically waits for A and B to fully finish before running.
#
# IMPORTANT: run the length calibration FIRST (see calibrate below) before
# running this - TOTAL_LENGTH_MM in build_design.py must already be
# corrected to hit 5 GHz, since this script uses that module constant for
# every task.
#
# Run this ONCE from the frontend, then walk away. Check back later with:
#   qstat
# and once everything shows as finished:
#   cat ~/circuit_design/batch_test/stageC.o<jobid>
#   ls ~/sweep_v2/plots/
#
# Usage:
#   ./run_full_sweep.sh            # full 10-position sweep
#   ./run_full_sweep.sh --pilot    # just tasks 1-2 (baseline + 1 numeric
#                                   # position) of Stage A/B, to validate
#                                   # the whole pipeline end-to-end cheaply
#                                   # before committing 10x the cluster
#                                   # time/memory risk to the full sweep.
#                                   # Stage C still waits for and
#                                   # aggregates whatever ran.

set -e
cd "$(dirname "$0")"

TASK_RANGE=""
if [[ "$1" == "--pilot" ]]; then
    echo "=== PILOT MODE: running only tasks 1-2 (baseline + 1 numeric position) ==="
    echo "    Use this to validate the full pipeline (mesh settings, node "
    echo "    targeting, wall-clock pace) before committing to all 10 "
    echo "    positions. Re-run without --pilot once satisfied."
    echo ""
    TASK_RANGE="-t 1-2"
fi

echo "=== Submitting Stage A (eigenmode sweep) ==="
qsub $TASK_RANGE stageA_eigenmode.sh

echo "=== Submitting Stage B (driven S21 sweep) ==="
qsub $TASK_RANGE stageB_driven.sh

echo "=== Submitting Stage C (aggregation + plots, waits for A and B) ==="
qsub stageC_aggregate.sh

echo ""
echo "=== All jobs submitted ==="
echo "Check progress with:  qstat"
echo "  - stageA / stageB tasks will show state 'qw' (queued) then 'r' (running)"
echo "  - stageC will show 'hqw' (held, queued) until ALL stageA+stageB tasks finish"
echo ""
echo "Check which node a running task landed on (important - see stageB's"
echo "own comments on node-size requirements):"
echo "  qstat -t     # shows queue@hostname per task"
echo ""
echo "When stageC has run and disappeared from qstat, check:"
echo "  cat $(pwd)/stageC.o*"
echo "  ls -la \$HOME/sweep_v2/plots/"
echo ""
echo "If anything failed, check individual task logs:"
echo "  cat $(pwd)/stageA.o<jobid>.<taskid>"
echo "  cat $(pwd)/stageB.o<jobid>.<taskid>"
echo ""
echo "IMPORTANT: check the ground_pos=none (baseline) task log first, since"
echo "it validates whether TOTAL_LENGTH_MM in build_design.py is actually"
echo "calibrated correctly for a 5 GHz fundamental - look for the line:"
echo '  *** BASELINE CHECK: measured f1=... vs target 5.0000 GHz ... ***'
echo "in stageA task 1 log."
