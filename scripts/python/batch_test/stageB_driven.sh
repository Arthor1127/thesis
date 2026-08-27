#!/bin/bash
#$ -cwd
#$ -N stageB
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=150G
#$ -l h=!compute-4-0&!compute-4-1&!compute-4-2&!compute-4-3&!compute-4-5&!compute-4-6&!compute-4-7&!compute-4-8&!compute-4-9&!compute-4-11&!compute-4-12&!compute-4-13
#$ -pe neworte 16
#$ -t 1-10
#$ -tc 2
#$ -V

# Each task runs ONE continuous 4.5-20.5 GHz sweep (covering all four
# harmonics at once) for one grounding position - no per-harmonic
# splitting anymore, so 10 positions = 10 tasks.
#
# Does not depend on Stage A (fixed frequency grid, not per-point
# resonance), so this can run concurrently with Stage A.
#
# --- Queue (copahue,highmem) ---
# compute-4-18 (503.4G, by far the largest node available - and the one
# most worth using given this stage's memory profile) is NOT a member of
# the copahue queue by itself - qhost -q confirmed it only appears under
# 'highmem' and 'be_copahue'. With -q copahue alone it's invisible to
# the scheduler regardless of the h= exclusions below, and was observed
# sitting fully idle for 14+ hours while every OTHER eligible node was
# full. Listing both queues here fixes that.
#
# --- Node targeting (memoria_a_usar=150G + explicit exclusions) ---
# A single confirmed production driven run (order=2, split mesh,
# ~42M unknowns) pushed a 125.7G node (compute-4-1) to MEMUSE=122.5G
# and rising swap (14.7G) - i.e. even a "normal-sized" node was NOT
# safely large enough. The exclusion list above removes every node
# under ~188G (compute-4-0/1/2/3/5/6/7/9/11/12/13, all <=141.4G),
# leaving only compute-4-4 (220.2G), compute-4-14/15/16/17 (188.4G
# each), and compute-4-18 (503.4G) eligible. Re-check `qhost` if this
# ever needs revisiting - node sizes/availability can change.
#
# --- Adaptive fast frequency sweep ---
# --adaptive-tol 1e-3 switches Palace from solving the full-order model
# at all 321 frequencies to solving it at a few auto-chosen samples
# (<=30) and evaluating the rest from a reduced-order model. Since each
# full-order solve rebuilds the AMG preconditioner and that is ~75% of
# runtime, this is potentially a ~10x cut in Stage B wall time. Validate
# the first position against a --adaptive-tol 0 run before trusting all
# ten: the tolerance is an L2 indicator, not a bound on S-parameters.
#
# --- Cores (16, NOT 32) ---
# An earlier note here claimed 32 cores was a measured 8x win. That was
# wrong: the two runs compared (4117s and 494s) differed in BOTH core
# count AND mesh (3.91M vs 0.33M unknowns), so they measure the mesh, not
# the cores. Worse, time improved only 8.3x while DOF dropped 11.8x, so
# there is no evidence the extra cores helped at all.
#
# Scheduling also argues against 32. Only compute-4-14/15/16/17 (copahue)
# and compute-4-18 (highmem) offer 32 slots, and the copahue four were
# observed persistently full with other users' multi-week jobs - a -pe 32
# job realistically competes for compute-4-18 alone, so -tc 2 would still
# only ever place one task. At 16 slots both tasks can run, and more
# nodes qualify. Throughput beats per-job latency when there are 10
# positions to get through.
#
# If you want to settle the core-count question, run the SAME mesh at 16
# and 32 cores and compare the Preconditioner line in the timing report.
#
# --- Concurrency throttle (-tc 2) ---
# The array has 10 tasks but no per-node memory isolation is
# guaranteed by SGE here (memoria_a_usar has been observed to NOT be
# strictly enforced - an earlier job landed on a too-small node
# despite requesting more than it had). -tc 2 caps how many stageB
# tasks run at once CLUSTER-WIDE, reducing the chance of multiple
# ~100-125G+ jobs landing on the same large node simultaneously and
# collectively exceeding it. This trades wall-clock time (a full
# sweep will take longer than max parallel) for memory safety - raise
# this only after confirming actual per-run peak memory on THIS mesh/
# solver-order combination via `qacct -j <jobid>` on a completed task.

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

POSITIONS=(none $(python3 -c "from build_design import sweep_positions_mm; print(' '.join(str(x) for x in sweep_positions_mm(n_points=9)))"))
# 1 baseline + 9 numeric positions = 10 tasks -> matches #$ -t 1-10 above.

POS=${POSITIONS[$((SGE_TASK_ID-1))]}

OUTDIR=$HOME/sweep_v2/driven_${POS}
mkdir -p $OUTDIR

echo "=== Stage B task $SGE_TASK_ID: ground_pos=$POS mm, continuous 4.5-20.5 GHz sweep ==="
echo "=== Running on host: $(hostname) ==="
python build_and_driven.py \
    --ground-pos $POS \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --bg-min-size-um 15 --bg-max-size-um 150 \
    --resonator-bg-min-size-um 30 --resonator-bg-max-size-um 300 \
    --solver-order 2 \
    --solns-to-save 10 \
    --solver-tol 1e-8 \
    --adaptive-tol 1e-3 \
    --adaptive-max-samples 30
