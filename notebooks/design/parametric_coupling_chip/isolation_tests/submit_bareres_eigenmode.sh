#!/bin/bash
#$ -N qc_bareres
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue,highmem,be_copahue
#$ -l memoria_a_usar=24G
#$ -pe neworte 16
#$ -o bareres_task_$TASK_ID.log
#$ -t 1-6

# Isolation test: bare shorted quarter-wave CPW resonator, ONE lumped port
# (resistance R) at the open end, no coupler, no launchpads -- an analytic
# Qc = (pi/4)*(R/Z0) exists to check the eigenmode-lossy-port method against
# with zero coupler-geometry ambiguity. See run_bareres_eigenmode_isolation.py
# for the full motivation.
#
# Default PARAMS_FILE (params_bareres_eigenmode.txt) is the original 3-point
# Qc~3k/16k/79k bracket. For the ceiling investigation pass
#   -v SOLVER_TOL=1e-8,PARAMS_FILE=isolation_tests/params_bareres_eigenmode_tol8.txt
# (6 rows, Qc_analytic ~82k..1e6). At SOLVER_TOL=1e-8/order 1 the hard 4.5e7
# saturation cliff is gone but loaded_Q under-reads by ~15% at Qc=150k rising
# to ~70% at Qc=1e6. The order-4 pass targets that residual bias:
#   qsub -q highmem -tc 2 -l memoria_a_usar=180G \
#        -v SOLVER_TOL=1e-8,SOLVER_ORDER=4,SOLVER_MAXITS=300,NUM_FREQS=4,\
#           PARAMS_FILE=isolation_tests/params_bareres_eigenmode_tol8.txt \
#        isolation_tests/submit_bareres_eigenmode.sh
# (order 4 is ~10-50x the RAM/time of order 1 -- needs the 503G highmem node).
# Solver knobs: SOLVER_TOL sets BOTH Eigenmode.Tol and Linear.Tol in SQDMetal;
# SOLVER_ORDER sets Solver.Order; SOLVER_MAXITS sets Linear.MaxIts; NUM_FREQS
# sets Eigenmode.N. All fall back to the driver's own defaults if not passed.

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

# Stagger array-task starts: two tasks hitting SQDMetal's setup in the shared
# cwd in the same second raced on the PALACE output folder once (task 2 of
# job 1938972 died with "Error creating test file in output folder").
sleep $(( (SGE_TASK_ID - 1) * 8 ))

PARAMS_FILE="${PARAMS_FILE:-isolation_tests/params_bareres_eigenmode.txt}"
LINE=$(sed -n "${SGE_TASK_ID}p" "$PARAMS_FILE")
read -r R TAG <<< "$LINE"
export PORT_R_OHM="$R" TAG PALACE_BIN
export STARTING_FREQ_GHZ=7.5
export NUM_FREQS="${NUM_FREQS:-2}"
export SOLVER_TOL="${SOLVER_TOL:-1e-6}" SOLVER_MAXITS="${SOLVER_MAXITS:-100}" SOLVER_ORDER="${SOLVER_ORDER:-1}"
export LINEAR_TYPE="${LINEAR_TYPE:-}" KSP_TYPE="${KSP_TYPE:-}"
export OUTDIR="${OUTDIR:-out_bareres_eig}"
mkdir -p "$OUTDIR"

# RAM poller (ps --ppid only sees the shell's direct children, NOT the
# MPI-spawned palace ranks -- the trustworthy number is palace's own
# "Estimated ... memory ..." stdout line, grepped into the exit echo below).
( while true; do printf "%s %s\n" "$(date +%s)" \
    "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
    sleep 5; done ) > "$OUTDIR/ram_${SGE_TASK_ID}.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

RUN_LOG="$OUTDIR/run_${SGE_TASK_ID}.log"
echo "=== qc_bareres task $SGE_TASK_ID : R=$R -> TAG=$TAG host=$(hostname) NSLOTS=$NSLOTS solver_tol=$SOLVER_TOL solver_order=$SOLVER_ORDER solver_maxits=$SOLVER_MAXITS num_freqs=$NUM_FREQS $(date) ==="
python -u isolation_tests/run_bareres_eigenmode_isolation.py 2>&1 | tee "$RUN_LOG"
ST=${PIPESTATUS[0]}
kill $MON 2>/dev/null
echo "=== qc_bareres task $SGE_TASK_ID ($TAG) exit=$ST  peak RSS(shell): $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' "$OUTDIR/ram_${SGE_TASK_ID}.log" 2>/dev/null)  palace: $(grep -iE 'estimated.*memory|peak.*memory' "$RUN_LOG" | tail -1)  $(date) ==="
exit $ST
