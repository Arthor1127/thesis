#!/bin/bash
#$ -N qc_meshbr
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=64G
#$ -pe neworte 16
#$ -o mesh_bracket_$TASK_ID.log
#$ -t 1-2

# AMR does not work for Driven sims in this PALACE build (confirmed: jobs
# 1938861 and 1938862 both segfault at rank 0, ~70-90s into the solve,
# identically with and without PROM -- see run_qc_sim.py's AMR block and this
# session's own investigation). Falling back to filling in the manual mesh
# bracket between baseline (mesh_factor=1.0, Qc=14078) and mesh_finer
# (mesh_factor=0.5, Qc=15744, +11.8%) with intermediate points, to see
# whether Qc is converging smoothly (extrapolatable) or erratically as the
# mesh refines toward the coupling gap's true resolved field.
#
# params_meshbracket.txt: (MESH_FACTOR, TAG). 64G requested as a safety
# margin above baseline's ~15-25G, well under mesh_finer's surprising 103.6G
# -- if either point needs more, that's itself useful information about how
# fast memory is actually scaling with mesh_factor.

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
mkdir -p out_meshbr

LINE=$(sed -n "${SGE_TASK_ID}p" convergence_studies/params_meshbracket.txt)
read -r MESH_FACTOR TAG <<< "$LINE"
export COUPLING_SPACE_UM=14 COUPLING_LENGTH_UM=300
export FSTART_GHZ=7.6 FEND_GHZ=8.35 NPTS=12001
export ADAPTIVE_MAX_SAMPLES=300 SOLVER_TOL=1e-6 SOLVER_MAXITS=500
export MESH_FACTOR TAG PALACE_BIN OUTDIR="out_meshbr"

echo "=== qc_meshbr task $SGE_TASK_ID ($TAG): mesh_factor=$MESH_FACTOR host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u driven_sweep/run_qc_sim.py 2>&1
ST=$?
echo "=== qc_meshbr task $SGE_TASK_ID ($TAG) exit=$ST $(date) ==="
exit $ST
