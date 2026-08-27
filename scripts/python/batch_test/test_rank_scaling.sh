#!/bin/bash
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q highmem
#$ -l h=compute-4-18
#$ -l memoria_a_usar=64G
#$ -pe neworte 16
#$ -V

# RANK-SCALING EXPERIMENT
#
# Measures how Palace actually scales with MPI rank count on THIS
# problem, so the choice of -pe is based on data instead of assumption.
# Earlier in this project a "32 cores is 8x faster" claim turned out to
# be a mesh-size effect, because the two runs compared differed in mesh
# AND rank count. This script exists to not repeat that.
#
# HOW TO RUN - submit the SAME script three times, overriding -pe on the
# command line (a command-line -pe overrides the directive above):
#
#     qsub -pe neworte 8  -N rank08 test_rank_scaling.sh
#     qsub -pe neworte 16 -N rank16 test_rank_scaling.sh
#     qsub -pe neworte 32 -N rank32 test_rank_scaling.sh
#
# then summarise with:
#
#     python parse_rank_scaling.py rank08.o* rank16.o* rank32.o*
#
# WHY IT IS PINNED TO ONE HOST (-q highmem -l h=compute-4-18):
# comparing runs that land on different nodes measures the hardware, not
# the rank count. compute-4-18 is the only node reachable with 32 slots
# that is reliably idle (503G, 32 cores, and NOT a member of copahue -
# which is why the earlier stageB jobs queued 14h without ever seeing
# it). All three runs therefore execute on identical silicon.
#
# EVERYTHING ELSE IS HELD FIXED: same design, same length, same 8/100um
# mesh, same solver order, same shift target. Only $NSLOTS varies.
#
# WHAT TO LOOK FOR in the Elapsed Time Report:
#   * Preconditioner - the line that matters. It is ~73% of runtime, so
#     it alone sets the Amdahl ceiling: even made infinitely fast, total
#     speedup caps at ~3.7x.
#   * Total - the honest end-to-end number.
#   * If Preconditioner scales poorly from 16->32, that is evidence the
#     ghost-cell duplication and multigrid-coarsening fragmentation from
#     more subdomains are eating the gain. That would ALSO be the
#     argument for a hybrid MPI+OpenMP build (PALACE_WITH_OPENMP), which
#     cuts subdomain count for the same core count. If it scales well,
#     the hybrid build is not worth the effort.
#
# Cost: ~510s per run at 16 ranks, so the whole experiment is well under
# an hour of wall time.

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

# Outdir keyed on rank count so the three runs cannot overwrite each other.
OUTDIR=$HOME/sweep_qw/rankscale_${NSLOTS}
mkdir -p $OUTDIR

echo "=== RANK SCALING: NSLOTS=$NSLOTS ==="
echo "=== Host: $(hostname)  (pinned - all runs must show compute-4-18) ==="
echo "=== Design: quarter-wave baseline, L=5.92905mm, mesh 8/100um ==="

python build_and_eigenmode_qw.py \
    --ground-pos none \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --total-length-mm 5.92905 \
    --starting-freq 5.0e9 \
    --bg-min-size-um 8 --bg-max-size-um 100 \
    --target-freq-ghz 5.0
