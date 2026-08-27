#!/bin/bash
#$ -cwd
#$ -N qw_sweep
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=16G
#$ -pe neworte 16
#$ -t 1-10
#$ -tc 6
#$ -V

# End-to-end grounding-strap sweep on the BARE quarter-wave resonator
# (no TL, no launchpads). Task 1 is the ungrounded baseline; tasks 2-10
# walk the strap from x0=0.5mm to x0=5.429mm, measured from the SHORT
# end (sg1), since that is where RouteMeanderGrounded measures
# 'positions' from - it uses the route's start_pin, and
# build_design_quarterwave.py wires start_pin=sg1.
#
# --- Why each task gets its OWN --starting-freq ---
# A strap at x0 splits the line into two independent resonators:
#     A: [0,x0]  short-short -> n*v/(2*x0)
#     B: [x0,L]  short-open  -> (2m-1)*v/(4*(L-x0))
# so the lowest mode moves across the sweep, and not monotonically:
# 5.46 GHz near the short end, peaking ~14.1 GHz around x0/L=0.71, back
# to 10.9 GHz near the open end. Palace's eigensolver is shift-and-
# invert, so a single fixed target would sit far from the true mode at
# most positions. On the half-wave design that exact mistake produced an
# eig.csv with a header and ZERO converged modes. qw_sweep_plan.py
# predicts each task's target analytically.
#
# The predictions assume an IDEAL hard short. A real 10um strap has
# finite width and inductance, so these are search CENTRES, not expected
# answers - real modes should land somewhat below. --number-of-freqs 6
# gives the solver room to return neighbours either side, and
# build_and_eigenmode_qw.py's --target-freq-ghz then picks the closest
# one rather than trusting a fixed list index (a strap restructures the
# mode ladder, so index-based selection is not stable across a sweep).
#
# --- Point of interest ---
# x0 = 2L/3 = 3.9527mm sits between tasks 7 (3.581) and 8 (4.197). There
# BOTH segments resonate at 3*f1 = 15 GHz simultaneously - a degenerate
# pair - because that is the 3*lambda/4 mode's own node, so the strap
# barely perturbs it. This coarse sweep will straddle but not land on
# it; run a dense sweep around 3.9527mm afterwards to resolve the
# degeneracy properly.
#
# --- Resources, re-derived after the fine_mesh_in_rectangle max_size fix ---
# Everything below was originally sized for meshes inflated ~36x by that
# bug (the rectangle field's max_size was set to the background MIN, which
# capped the WHOLE domain at that size and destroyed far-field coarsening;
# a strapped mesh measured 7.81M elements at EVERY position, vs 214,924
# after the fix). Peak memory on this design has since run ~10G, so:
#   memoria_a_usar 32G -> 16G
#   node exclusion list dropped entirely - at this footprint even the
#     small 30-31G nodes qualify, and restricting to a handful of hosts
#     was what left an earlier job queued 14h while compute-4-18 sat idle
#   -tc 4 -> 6, since each task is far cheaper and 6*16=96 slots is well
#     inside the 200-slot quota confirmed by qquota
#
# --- Settings ---
# 16 cores (not 32): the 32-core "win" measured earlier on the half-wave
# design was a mesh-size artifact, not a core-count effect - see
# stageA_eigenmode.sh. 8/100um mesh: the measured accuracy/cost sweet
# spot; never go above ~6um min or the mesh stops resolving the CPW
# cross-section (a 20um run produced spurious clustered modes). No AMR:
# --amr-max-its exists but is unvalidated, do not enable blind.

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg

PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/batch_test"
cd "$SCRIPT_DIR"

# Params for THIS task, from the two-segment model. Parsed without
# process substitution so this does not depend on the shell supporting
# it, and guarded so a planner failure aborts the task rather than
# silently running build_and_eigenmode_qw.py with empty arguments.
PARAMS=$(python qw_sweep_plan.py --task-id $SGE_TASK_ID) || {
    echo "FATAL: qw_sweep_plan.py failed for task $SGE_TASK_ID" >&2
    exit 1
}
POS=$(echo "$PARAMS" | awk '{print $1}')
START_FREQ=$(echo "$PARAMS" | awk '{print $2}')
TARGET_GHZ=$(echo "$PARAMS" | awk '{print $3}')

if [ -z "$POS" ] || [ -z "$START_FREQ" ] || [ -z "$TARGET_GHZ" ]; then
    echo "FATAL: could not parse task params, got: '$PARAMS'" >&2
    exit 1
fi

OUTDIR=$HOME/sweep_qw/geom_${POS}
mkdir -p $OUTDIR

echo "=== QW sweep task $SGE_TASK_ID: ground_pos=$POS mm (from SHORT end) ==="
echo "=== shift-invert target = $TARGET_GHZ GHz (predicted lowest mode) ==="
echo "=== Running on host: $(hostname) ==="

python build_and_eigenmode_qw.py \
    --ground-pos $POS \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --total-length-mm 5.92905 \
    --starting-freq $START_FREQ \
    --number-of-freqs 6 \
    --bg-min-size-um 8 --bg-max-size-um 100 \
    --target-freq-ghz $TARGET_GHZ
