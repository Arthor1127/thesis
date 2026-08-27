#!/bin/bash
#$ -cwd
#$ -N eig_fast
#$ -j y
#$ -S /bin/bash
#$ -q highmem
#$ -l memoria_a_usar=100G
#$ -pe neworte 32
#$ -V

# FAST calibration run - same question as test_eig_finemesh.sh (what
# TOTAL_LENGTH_MM gives a true 5GHz fundamental?) but aiming for a
# fraction of that run's 4117s / 69min wall time.
#
# What changed, and why, based on the 6/80um run's own timing report:
#
#   Preconditioner   3002s   <-- 75% of total runtime
#   Linear Solve      349s
#   Coarse Solve      235s
#   Div-Free Proj     167s
#   Estimation Solve  137s
#   everything else  ~120s
#
# Runtime is dominated by AMG preconditioner SETUP, which scales worse
# than linearly with unknown count. So the lever is fewer unknowns, not
# a faster solve.
#
# 1) -pe neworte 32 (was 16). compute-4-18 has 32 cores and the highmem
#    queue grants all 32; the 6/80um run peaked at 32.2G of 503G, so
#    there is enormous headroom. AMG setup parallelises reasonably well.
#    Free ~1.5-2x. Queue is highmem-only here since that is the only
#    queue offering 32 slots on this node.
#
# 2) --mesh-refinement 3 (was 0, i.e. AMR off). Palace was ALREADY
#    computing the error indicator every run and then discarding it.
#    With AMR on it refines only where the indicator is large, instead
#    of the uniform global refinement the 15->8->6um sweep was doing.
#    Same accuracy at far fewer DOF.
#
# 3) Coarse STARTING mesh (20/200um). Deliberately coarser than any
#    previous run - AMR adds resolution where the physics needs it, so
#    starting fine just wastes elements in the empty substrate.
#
# 4) --starting-freq 5.0e9, since we now know f1 sits within ~1% of
#    5GHz at this length. A tight shift-invert target converges faster
#    than a distant one.
#
# EXPECTATION, not a guarantee: this should land at a similar or better
# f1 than the 6/80um run's 5.0389 GHz in materially less wall time. If
# AMR turns out to add more setup cost than it saves, that shows up as
# a longer runtime and is itself a useful result - compare the elapsed
# time report against the numbers in the table above.

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

OUTDIR=$HOME/sweep_v2/eig_fast_test
mkdir -p $OUTDIR

echo "=== FAST eigenmode calibration: ground_pos=none, L=11.7914mm ==="
echo "=== 32 cores, AMR=3, coarse 20/200um start mesh ==="
echo "=== Running on host: $(hostname) ==="

python build_and_eigenmode.py \
    --ground-pos none \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --total-length-mm 11.7914 \
    --starting-freq 5.0e9 \
    --bg-min-size-um 20 --bg-max-size-um 200 \
    --mesh-refinement 3 \
    --target-freq-ghz 5.0
