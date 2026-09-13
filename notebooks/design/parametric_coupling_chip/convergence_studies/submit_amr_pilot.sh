#!/bin/bash
#$ -N qc_amr
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=48G
#$ -pe neworte 16
#$ -o amr_pilot.log

# AMR pilot at the cs14_cl300 reference point (fr=8.101GHz, Qc=14078 from the
# manual, non-AMR baseline). Base mesh forced to mesh_factor=1.0 (see
# run_qc_sim.py's own enforcement -- AMR cannot rescue a too-coarse starting
# mesh, only sharpen an adequate one). 3 iterations, DOF capped at 1.5M as a
# first conservative guess (mesh_finer, a single uniform 2x refinement of
# this same base mesh, needed 103.6GB -- watch this pilot's own peak-memory
# report closely and do not scale AMR_MAX_DOFS up without evidence it's safe).
# 16 cores (not 8): AMR reruns the solve per iteration, so wanted more headroom
# than a single-solve convergence-sweep point; still well under a full-node
# request.

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
mkdir -p out_amr

export COUPLING_SPACE_UM=14 COUPLING_LENGTH_UM=300
export FSTART_GHZ=7.6 FEND_GHZ=8.35 NPTS=12001
export ADAPTIVE_MAX_SAMPLES=300 SOLVER_TOL=1e-6 SOLVER_MAXITS=500
export AMR_ITERS=3 AMR_TOL=1e-2 AMR_MAX_DOFS=1500000 AMR_DORFLER=0.7
export TAG=amr_pilot PALACE_BIN OUTDIR="out_amr"

echo "=== qc_amr pilot host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u driven_sweep/run_qc_sim.py 2>&1
ST=$?
echo "=== qc_amr pilot exit=$ST $(date) ==="
exit $ST
