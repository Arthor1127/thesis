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
#   ./run_full_sweep.sh

set -e
cd "$(dirname "$0")"

echo "=== Submitting Stage A (eigenmode sweep, 10 tasks) ==="
qsub stageA_eigenmode.sh

echo "=== Submitting Stage B (driven S21 sweep, 10 tasks) ==="
qsub stageB_driven.sh

echo "=== Submitting Stage C (aggregation + plots, waits for A and B) ==="
qsub stageC_aggregate.sh

echo ""
echo "=== All jobs submitted ==="
echo "Check progress with:  qstat"
echo "  - stageA / stageB tasks will show state 'qw' (queued) then 'r' (running)"
echo "  - stageC will show 'hqw' (held, queued) until ALL stageA+stageB tasks finish"
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
