#!/bin/bash
#$ -N qc_conv
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=28G
#$ -pe neworte 8
#$ -tc 6
#$ -o conv_task_$TASK_ID.log
#$ -t 1-14

# RELIABILITY (not speed) convergence study: does run_qc_sim.py's fitted
# (fr, Qc) actually agree/stabilize as solver_tol, solver_maxits, mesh
# resolution, and AdaptiveMaxSamples are bracketed both looser AND tighter
# than our current defaults? If the answer drifts under refinement, the
# current numbers can't be trusted regardless of runtime; if it's stable
# across the bracket, that's the evidence needed to trust them going forward.
# See params_convergence.txt for the (SOLVER_TOL, SOLVER_MAXITS, MESH_FACTOR,
# ADAPTIVE_MAX_SAMPLES, TAG) grid -- one parameter varied at a time from the
# 'baseline' row (row 1), which reproduces run_qc_sim.py's exact current
# defaults exactly (FSTART/FEND/NPTS unchanged too).
#
# Fixed point: coupling_space_um=14, coupling_length_um=300 (the same point
# already carefully characterized this session: fr=8.101GHz, Qc=14078,
# Qr=533). Window deliberately kept at the full, already-validated
# 7.6-8.35GHz/12001-point setting (NOT narrowed to save cost) -- this sweep
# runs unattended overnight, and a narrower window risks missing the
# resonance entirely if a looser solver/coarser mesh setting shifts the
# apparent fr enough to fall outside it. A missed resonance under a
# DELIBERATELY-relaxed setting is exactly the failure mode this reliability
# check needs to be able to see and report, not accidentally mask by having
# already cropped the search window down to where we know the answer sits.
#
# Queue: copahue only (NOT highmem -- a parallel agent uses that queue for the
# same kind of jobs, per standing instruction). 8 cores/task, small-core
# fan-out across many points, per this repo's own resource-sizing convention.

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
mkdir -p out_conv

PARAMS_FILE="${PARAMS_FILE:-convergence_studies/params_convergence.txt}"
LINE=$(sed -n "${SGE_TASK_ID}p" "$PARAMS_FILE")
read -r SOLVER_TOL SOLVER_MAXITS MESH_FACTOR ADAPTIVE_MAX_SAMPLES TAG <<< "$LINE"
export COUPLING_SPACE_UM=14 COUPLING_LENGTH_UM=300
export FSTART_GHZ=7.6 FEND_GHZ=8.35 NPTS=12001
export SOLVER_TOL SOLVER_MAXITS MESH_FACTOR ADAPTIVE_MAX_SAMPLES TAG PALACE_BIN
export OUTDIR="out_conv"

echo "=== qc_conv task $SGE_TASK_ID ($TAG): tol=$SOLVER_TOL maxits=$SOLVER_MAXITS mesh_factor=$MESH_FACTOR ams=$ADAPTIVE_MAX_SAMPLES host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u driven_sweep/run_qc_sim.py 2>&1
ST=$?
echo "=== qc_conv task $SGE_TASK_ID ($TAG) exit=$ST $(date) ==="
exit $ST
