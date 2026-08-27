#!/bin/bash
#$ -cwd
#$ -N stageA
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=64G
#$ -l h=!compute-4-6&!compute-4-7&!compute-4-8
#$ -pe neworte 16
#$ -t 1-10
#$ -tc 4
#$ -V

# --- Queue (copahue,highmem) ---
# compute-4-18 (503.4G, by far the largest node available) is NOT a
# member of the copahue queue - qhost -q confirmed it only appears under
# 'highmem' and 'be_copahue'. With -q copahue alone, this node is
# invisible to the scheduler regardless of h= exclusions, and was
# observed sitting completely idle while stageB queued for 14+ hours
# because its other 5 "eligible" nodes were all full. Listing both
# queues here lets SGE actually schedule onto compute-4-18 when needed.

# --- Slot quota (-tc 4) ---
# Per the cluster docs (instrucciones-utilizar-cluster.html), personal
# quota is 140 slots combined across caulle+copahue queues. Without a
# throttle here, this array alone could request up to 10*16=160 slots -
# already over quota by itself, and would compete with stageB_driven.sh
# (capped at -tc 2 = 32 slots) for the same personal allocation if both
# run concurrently (run_full_sweep.sh submits both independently). -tc 4
# caps this at 4*16=64 slots. Paired with stageB's -tc 2 * 16 = 32, that
# is 96 of the 140 limit with both stages running - comfortable.


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

# Positions generated from build_design.py's own sweep_positions_mm() so
# this never drifts out of sync with the design module. "none" (the
# ungrounded baseline) is always task 1.
POSITIONS=(none $(python3 -c "from build_design import sweep_positions_mm; print(' '.join(str(x) for x in sweep_positions_mm(n_points=9)))"))
# 1 baseline + 9 numeric positions = 10 tasks -> matches #$ -t 1-10 above.
# If you change n_points in sweep_positions_mm(), update the -t range to
# (n_points + 1).

POS=${POSITIONS[$((SGE_TASK_ID-1))]}

OUTDIR=$HOME/sweep_v2/geom_${POS}
mkdir -p $OUTDIR

echo "=== Stage A task $SGE_TASK_ID: ground_pos=$POS mm ==="
echo "=== Running on host: $(hostname) ==="
# Mesh floor: --bg-min-size-um must stay at or above ~6um. The CPW has
# a 10um trace and 6um gaps; a 20um run produced f1=6.11GHz and clustered
# spurious modes because the mesh could not resolve the cross-section at
# all. 8/100um is the accuracy/cost sweet spot measured so far.
#
# TOTAL_LENGTH_MM in build_design.py is still pending a final decision -
# see LESSONS_LEARNED. Normalising f*L across the 15/150, 8/100 and 6/80
# runs suggested the mesh-converged length sits near ~11.97mm, close to
# the original elliptic-integral prediction, so some of the length
# iteration was chasing mesh error rather than physics.
python build_and_eigenmode.py \
    --ground-pos $POS \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --bg-min-size-um 8 --bg-max-size-um 100 \
    --target-freq-ghz 5.0
