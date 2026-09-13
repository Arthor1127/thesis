#!/bin/bash
#$ -N qc_lineConv
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem,copahue
#$ -l memoria_a_usar=32G
#$ -pe neworte 8
#$ -tc 4
#$ -o lineConv_task_$TASK_ID.log
#$ -t 1-7

# Mesh-resolution x vacuum-height (z_prop_top) convergence sweep on the
# line_only bare-CPW-line geometry, to find out why the simulated Zin (~33
# Ohm, see out_line/results/line_only.json) disagrees with BOTH DesignHelper's
# and SQDMetal's own CpwParams analytic formula (both ~50 Ohm). Traced the
# domain construction in SQDMetal.PALACE.Model.PALACE_Model_Base.__init__:
# defaults are set_xBoundary_as_proportion(0.1), set_yBoundary_as_proportion
# (0.1), set_zBoundary_as_proportion(1.0, 0.0) -- i.e. only 500um of vacuum
# above the chip (z_prop_top=1.0 x the 500um substrate) and a PEC floor
# exactly at the substrate's own bottom (z_prop_bottom=0, no margin).
# CpwParams.calc_impedance (a grounded-backplane-aware formula) predicts the
# bottom PEC at 500um is already converged (50.29 vs 50.31 Ohm at h->infinity)
# -- so that's not the likely culprit -- but it says nothing about whether
# 500um of vacuum ABOVE is enough, which is untested. This sweep checks BOTH
# vacuum height and mesh resolution, since an under-resolved mesh could
# independently cause an error of this size.
#
# See params_line_convergence.txt for the (MIN_SIZE_FACTOR, Z_PROP_TOP, TAG)
# grid: a mesh-only refinement series at the default z_prop_top=1.0, a
# z_prop_top-only series at the default mesh, and two combined-refinement
# points in case neither alone converges. All at MODE=line_only (cheapest,
# cleanest case -- no launchpad/coupler complexity).
#
# 8 cores x 4 concurrent: each of these ran in ~70-240s wall time at the
# default resolution: cheap, and this repo's own resource-sizing rule
# (cluster-resource-request-sizing memory) says fan many small-core jobs
# across the quota for exactly this kind of many-point-but-individually-light
# sweep rather than requesting large per-job core counts.

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
export OUTDIR="out_line_conv"
mkdir -p "$OUTDIR"

LINE=$(sed -n "${SGE_TASK_ID}p" convergence_studies/params_line_convergence.txt)
read -r MIN_SIZE_FACTOR Z_PROP_TOP TAG <<< "$LINE"
export MODE="line_only" MIN_SIZE_FACTOR Z_PROP_TOP TAG

echo "=== qc_lineConv task $SGE_TASK_ID ($TAG): min_size_factor=$MIN_SIZE_FACTOR z_prop_top=$Z_PROP_TOP host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u driven_sweep/check_line_impedance.py 2>&1
ST=$?
echo "=== qc_lineConv task $SGE_TASK_ID ($TAG) exit=$ST $(date) ==="
exit $ST
