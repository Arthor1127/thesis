#!/bin/bash
#$ -N qc_lineLat
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem,copahue
#$ -l memoria_a_usar=32G
#$ -pe neworte 8
#$ -tc 4
#$ -o lineLat_task_$TASK_ID.log
#$ -t 1-5

# Lateral ground-plane-width sweep, following up on submit_line_convergence.sh:
# that sweep found |Zin| climbing 33.3->38.0->40.5 Ohm as mesh refined
# (min_size_factor 0.5->0.25->0.125), but then PLATEAUING exactly at 40.483
# Ohm even under a further 2x mesh refinement AND a further 2x vacuum-height
# increase (conv_mesh4x_z6x vs conv_mesh8x_z12x -- identical to 3 decimals).
# So neither mesh resolution nor vacuum height explains the remaining
# 40.5-vs-50 Ohm gap. The one geometric truncation not yet tested: lateral
# ground-plane width. Defaults (x_prop=y_prop=0.1 on a 4.8x2mm chip) give only
# 200um of ground beyond the chip edge in y (total half-width ~1.2mm from the
# trace) and 480um in x -- this sweep tests whether widening that moves |Zin|
# further toward 50 Ohm. All points fixed at min_size_factor=0.125 (already
# shown equivalent to 0.0625, so no need to pay for the finer/slower mesh) and
# z_prop_top=1.0 (already shown irrelevant).
#
# See params_line_lateral.txt for the (MIN_SIZE_FACTOR, Z_PROP_TOP, X_PROP,
# Y_PROP, TAG) grid.

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
export OUTDIR="out_line_lateral"
mkdir -p "$OUTDIR"

LINE=$(sed -n "${SGE_TASK_ID}p" convergence_studies/params_line_lateral.txt)
read -r MIN_SIZE_FACTOR Z_PROP_TOP X_PROP Y_PROP TAG <<< "$LINE"
export MODE="line_only" MIN_SIZE_FACTOR Z_PROP_TOP X_PROP Y_PROP TAG

echo "=== qc_lineLat task $SGE_TASK_ID ($TAG): x_prop=$X_PROP y_prop=$Y_PROP host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u driven_sweep/check_line_impedance.py 2>&1
ST=$?
echo "=== qc_lineLat task $SGE_TASK_ID ($TAG) exit=$ST $(date) ==="
exit $ST
