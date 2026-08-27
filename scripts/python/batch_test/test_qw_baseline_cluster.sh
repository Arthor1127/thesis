#!/bin/bash
#$ -cwd
#$ -N eig_qw_baseline
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=32G
#$ -l h=!compute-4-6&!compute-4-7&!compute-4-8
#$ -pe neworte 16
#$ -V

# Quarter-wave BASELINE test - cluster, fast settings.
#
# Same design as test_qw_baseline_local.sh, submitted here instead of
# run locally. This design has no TL/launchpads (bare resonator only),
# so it should be substantially cheaper than any half-wave eigenmode run
# in this project - memoria_a_usar is set to a conservative 32G, well
# under what any run this session has needed, and node exclusions are
# the light half-wave-eigenmode set (not the strict driven-only list),
# since there is no reason to expect this needs a huge node.
#
# 16 cores, not 32: an earlier attempt at 32 cores for the half-wave
# design was NOT a genuine speedup (see stageA_eigenmode.sh's own
# comments - the "8x win" was actually a mesh-size effect, not cores;
# time improved LESS than DOF dropped in that comparison). No reason to
# expect a smaller design behaves differently, so stick with 16.
#
# Mesh: 8/100um - the measured accuracy/cost sweet spot for the
# half-wave design's CPW (10um trace/6um gap). This design's trace is
# wider (20um/12.25um gap), so 8um is if anything MORE conservative
# here (finer relative to the larger feature size) - should be safe,
# but if you want to push speed further once the baseline is confirmed,
# 10/120um is a reasonable next mesh point to try given the larger trace.
#
# No AMR (--amr-max-its 0, the default) - SQDMetal's own
# 'mesh_refinement' key is UniformLevels (uniform subdivision, NOT AMR);
# the real Palace AMR keys are exposed via --amr-max-its/--amr-tol in
# build_and_eigenmode.py but are UNVALIDATED - do not enable blind on a
# new geometry.

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

OUTDIR=$HOME/sweep_qw/eig_baseline
mkdir -p $OUTDIR

echo "=== Quarter-wave BASELINE: ground_pos=none, L=5.92905mm ==="
echo "=== Running on host: $(hostname) ==="

python build_and_eigenmode_qw.py \
    --ground-pos none \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --total-length-mm 5.92905 \
    --starting-freq 5.0e9 \
    --bg-min-size-um 8 --bg-max-size-um 100 \
    --target-freq-ghz 5.0
