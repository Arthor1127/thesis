#!/bin/bash
#$ -cwd
#$ -N qpdk_launcher_resonator
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=48G
#$ -pe neworte 32
#$ -V

# Real cluster run of the GDS+stackup pipeline (packages/gds-stackup-pipeline):
# a genuine qpdk-provided geometry (launcher + coupled quarter-wave CPW
# resonator), stackup grounded in qpdk's own published Si/500um constants,
# submitted end-to-end (mesh build + Palace eigenmode solve) as one job, per
# the existing squid_array_coupling/submit_smoke.sh convention. 48G/32 slots
# matches that smoke-test sizing -- no reason to request more without
# evidence.
#
# Meshing: fine_mesh_features(max_feature_size=20e-6), not
# fine_mesh_in_rectangle() -- the latter meshes the ENTIRE interior of its
# rectangle at mesh_min uniformly (GMSH's distance field for a filled region
# is zero throughout its interior, not just near the metal). A first attempt
# sized to the whole circuit footprint hit the timeout below without
# converging; fine_mesh_features is geometry-aware (only the actual narrow
# CPW trace/gap gets fine-meshed) and produced a ~3x smaller mesh that
# finished in ~44 minutes. See packages/gds-stackup-pipeline/README.md's
# "Gotchas" section for this and other pitfalls hit getting this run working.

set -uo pipefail

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/gds_pipeline_deploy/notebooks/design/gds_pipeline"
cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/cluster_run_1"
mkdir -p "$OUTDIR"

echo "=== qpdk launcher+resonator: host=$(hostname), NSLOTS=$NSLOTS, $(date) ==="

# Idempotent -- installs/updates the deployed package into this job's env.
# Never run standalone via bare SSH on the frontend; only ever inside a job.
pip install -e "$SCRIPT_DIR/../../../packages/gds-stackup-pipeline" 2>&1 | tail -5

timeout 3600 python -u run_qpdk_launcher_resonator.py \
    --outdir "$OUTDIR" \
    --num-cpus "$NSLOTS" \
    --palace-bin "$PALACE_BIN" \
    --run \
    2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}
echo "=== qpdk launcher+resonator done, exit=$STATUS, $(date) ==="
exit $STATUS
