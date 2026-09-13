#!/bin/bash
#$ -N qc_eig
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue,highmem,be_copahue
#$ -l memoria_a_usar=40G
#$ -pe neworte 16
#$ -tc 4
#$ -o eig_task_$TASK_ID.log
#$ -t 1-1

# Eigenmode (lossy-port, complex-eigenfrequency) Qc extraction -- see
# run_qc_eigenmode.py's docstring for the full motivation and method. This
# validates the method against cs14_cl300, the same reference point the
# driven-sim recipe (mesh_factor=0.75/tol=1e-4/maxits=100/ams=150) was
# validated on: target is fr ~ 8.1975GHz, Qc ~ 15684 (from the driven S21
# fit). A single eigenmode solve per point -- no frequency sweep needed --
# so this is CHEAPER than the driven sim, not just more reliable.
#
# See params_qc_eigenmode.txt for the (COUPLING_SPACE_UM, COUPLING_LENGTH_UM,
# TAG) grid -- currently just the one validation point.

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

sleep $(( (SGE_TASK_ID - 1) * 8 ))   # stagger: shared-cwd SQDMetal setup races otherwise

PARAMS_FILE="${PARAMS_FILE:-eigenmode_sweep/params_qc_eigenmode.txt}"
LINE=$(sed -n "${SGE_TASK_ID}p" "$PARAMS_FILE")
read -r CS CL TAG <<< "$LINE"
export COUPLING_SPACE_UM="$CS" COUPLING_LENGTH_UM="$CL" TAG PALACE_BIN
export STARTING_FREQ_GHZ="${STARTING_FREQ_GHZ:-8.0}" NUM_FREQS="${NUM_FREQS:-3}"
export MESH_FACTOR="${MESH_FACTOR:-0.75}"
export SOLVER_TOL="${SOLVER_TOL:-1e-8}" SOLVER_MAXITS="${SOLVER_MAXITS:-100}" SOLVER_ORDER="${SOLVER_ORDER:-1}"
export LINEAR_TYPE="${LINEAR_TYPE:-}" KSP_TYPE="${KSP_TYPE:-}"
export OUTDIR="${OUTDIR:-out_qc_eig}"
mkdir -p "$OUTDIR"

RUN_LOG="$OUTDIR/run_${SGE_TASK_ID}.log"
echo "=== qc_eig task $SGE_TASK_ID : cs=$CS cl=$CL -> TAG=$TAG host=$(hostname) NSLOTS=$NSLOTS tol=$SOLVER_TOL order=$SOLVER_ORDER $(date) ==="
python -u eigenmode_sweep/run_qc_eigenmode.py 2>&1 | tee "$RUN_LOG"
ST=${PIPESTATUS[0]}
echo "=== qc_eig task $SGE_TASK_ID ($TAG) exit=$ST  palace: $(grep -iE 'estimated peak.*node' "$RUN_LOG" | tail -1)  $(date) ==="
exit $ST
