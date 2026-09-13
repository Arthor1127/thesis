#!/bin/bash
#$ -N qc_final
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=40G
#$ -pe neworte 16
#$ -o final_task_$TASK_ID.log
#$ -t 1-2

# Concluding runs for the reliability/convergence investigation. Two checks
# in one array:
#
# 1. (cs14_cl300_final) Does the validated "safe efficiency" combo -- loose
#    solver tolerance (1e-4, down from 1e-6), fewer solver iterations (100,
#    down from 500), fewer PROM adaptive samples (150, down from 300) -- still
#    reproduce the CORRECTED mesh_factor=0.75 answer (fr=8.1975GHz,
#    Qc=15684) when all combined together? Each relaxation was only
#    validated ONE AT A TIME from the mesh_factor=1.0 baseline so far; this
#    is the first test of the actual combined "fast AND reliable" recipe
#    intended for the real Qc sweep.
#
# 2. (cs10_cl400_final) Does the same mesh_factor=0.75 correction (~+11%
#    Qc vs the mesh_factor=1.0 baseline, found at cs14_cl300) generalize to a
#    DIFFERENT design point? cs10_cl400 was previously only characterized at
#    mesh_factor=1.0 (manual refit: fr=8.15620GHz, Qc=13398).
#
# See params_conclude.txt for the (COUPLING_SPACE_UM, COUPLING_LENGTH_UM,
# MESH_FACTOR, SOLVER_TOL, SOLVER_MAXITS, ADAPTIVE_MAX_SAMPLES, TAG) grid.
# Window kept at the full conservative 7.6-8.35GHz/12001-point setting (both
# points' resonances have been seen to sit safely inside it even as mesh
# shifts fr somewhat).

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
mkdir -p out_final

LINE=$(sed -n "${SGE_TASK_ID}p" convergence_studies/params_conclude.txt)
read -r CS CL MESH_FACTOR SOLVER_TOL SOLVER_MAXITS ADAPTIVE_MAX_SAMPLES TAG <<< "$LINE"
export COUPLING_SPACE_UM="$CS" COUPLING_LENGTH_UM="$CL"
export FSTART_GHZ=7.6 FEND_GHZ=8.35 NPTS=12001
export MESH_FACTOR SOLVER_TOL SOLVER_MAXITS ADAPTIVE_MAX_SAMPLES TAG PALACE_BIN
export OUTDIR="out_final"

echo "=== qc_final task $SGE_TASK_ID ($TAG): cs=$CS cl=$CL mesh_factor=$MESH_FACTOR tol=$SOLVER_TOL maxits=$SOLVER_MAXITS ams=$ADAPTIVE_MAX_SAMPLES host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u driven_sweep/run_qc_sim.py 2>&1
ST=$?
echo "=== qc_final task $SGE_TASK_ID ($TAG) exit=$ST $(date) ==="
exit $ST
