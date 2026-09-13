#!/bin/bash
#$ -N qc_prod
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue,highmem,be_copahue,be_cpuINgpu
#$ -l memoria_a_usar=40G
#$ -pe neworte 8
#$ -tc 8
#$ -o prod_task_$TASK_ID.log
#$ -t 1-35

# Full production Qc sweep -- the actual deliverable this whole session's
# reliability investigation was in service of. Uses the VALIDATED recipe from
# tonight's convergence study rather than the old default settings:
#
#   MESH_FACTOR=0.75  (element size 0.75x the "half smallest feature" rule).
#       mesh_factor=1.0 UNDERSHOOTS the converged Qc by ~11-12% (checked at
#       cs14_cl300: 14078 -> 15684 going 1.0->0.75, then only +0.4% more
#       going all the way to 0.5x at 4x the compute cost) -- 0.75 is the
#       cheapest mesh that's actually close to converged.
#   SOLVER_TOL=1e-4, SOLVER_MAXITS=100 (down from 1e-6/500): checked to agree
#       with the tight baseline to 8+ significant figures at cs14_cl300 --
#       purely a speed win, zero accuracy cost, UNTIL ~1e-3 which crashes
#       PALACE outright (do not loosen further).
#   ADAPTIVE_MAX_SAMPLES=150 (down from 300): also agreed with baseline to
#       8+ figures once the resonance's rough location is already known,
#       which it is for every point in this grid (all within the same
#       7.6-8.35GHz window every prior pilot in this design family has used).
#
# 8 cores/task (not 16-32): per this repo's own sweep-throughput convention
# (small-core jobs fanned across the quota beats fewer big-core jobs). Queue
# list explicitly includes highmem and the best-effort queues for this
# overnight batch (per explicit instruction overriding the usual
# copahue-only default) -- copahue had ~49 free slots and highmem ~24 when
# this was launched; be_copahue/be_cpuINgpu showed ~0 free at that instant
# (suspended/disabled while their underlying nodes' regular queues are busy)
# but are left in the list since that can change over the course of the
# night. tc=8 uses most of the currently-free capacity without starving
# other users, and the whole 35-point grid can't all run at once regardless
# of core count. Some points may not
# finish by morning; qstat -u arturo.ruiz | grep qc_prod shows progress.
#
# EXPECT the shallow-dip guard to trip for most/all points (dip_depth_db >
# -3.0 -> Qr/Qc/fr left null): every point characterized this session so far
# has needed a manual +-100MHz refit around its own dip_freq_ghz, same as
# cs14_cl300/cs10_cl400 did. This is normal, not a failure -- refit_convergence
# .py's approach (glob results, refit each around its own recorded dip) is
# the template; a similar batch-refit pass over out_qc_prod/results/ is the
# expected next step once this array finishes.
#
# Queue list includes highmem/best-effort per explicit override for this
# overnight batch (see above) -- copahue-only is the usual standing default.
#
# PARAMS_FILE lets this script run a DIFFERENT points file without touching
# params_qc.txt while its own array still has queued tasks -- those read
# their row by SGE_TASK_ID whenever they actually start, potentially hours
# from now, so editing params_qc.txt mid-array would corrupt it (same
# mistake already made and fixed earlier this session).

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_chip"
cd "$SCRIPT_DIR"
mkdir -p out_qc_prod

PARAMS_FILE="${PARAMS_FILE:-driven_sweep/params_qc.txt}"
LINE=$(sed -n "${SGE_TASK_ID}p" "$PARAMS_FILE")
read -r CS CL TAG <<< "$LINE"
export COUPLING_SPACE_UM="$CS" COUPLING_LENGTH_UM="$CL" TAG PALACE_BIN
export FSTART_GHZ=7.6 FEND_GHZ=8.35 NPTS=12001
export MESH_FACTOR=0.75 SOLVER_TOL=1e-4 SOLVER_MAXITS=100 ADAPTIVE_MAX_SAMPLES=150
export OUTDIR="out_qc_prod"

echo "=== qc_prod task $SGE_TASK_ID : cs=$CS cl=$CL -> TAG=$TAG host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u driven_sweep/run_qc_sim.py 2>&1
ST=$?
echo "=== qc_prod task $SGE_TASK_ID ($TAG) exit=$ST $(date) ==="
exit $ST
