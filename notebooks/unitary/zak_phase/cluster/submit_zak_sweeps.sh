#!/bin/bash
# ===========================================================================
#  Submit a Zak-phase sweep array + a dependent plotting job.
#  Run ON the cluster frontend, from  ~/unitary/zak_phase/
#
#  Usage:  ./submit_zak_sweeps.sh <grid> <outdir> ["queue1,queue2,..."]
#    e.g.  ./submit_zak_sweeps.sh grid_omegar.json results_omegar
#
#  1-core tasks, ~15 min each, checkpointed (resumable). -tc caps concurrency.
#  The plot job holds for THIS array, then (re)plots every results_* dir it finds.
# ===========================================================================
set -e
GRID=${1:-grid_gammaomegar.json}
OUT=${2:-results_gammaomegar}
QUEUES=${3:-"cpuINgpu,copahue,be_cpuINgpu,be_copahue"}     # no highmem (other agent's)
TC=64
CONDA=/share/apps/miniconda3/etc/profile.d/conda.sh
ENV=qcg-quantum-design
HERE=$(pwd)
source ${CONDA}; conda activate ${ENV}
NSWEEPS=$(python -c "import json;print(len(json.load(open('${GRID}'))['sweeps']))")
mkdir -p logs

cat > script-array-${OUT} <<EOF
#!/bin/bash
#\$ -clear
#\$ -S /bin/bash
#\$ -N zak_${OUT}
#\$ -cwd
#\$ -t 1-${NSWEEPS}
#\$ -tc ${TC}
#\$ -l mem=1G
#\$ -j y
#\$ -o logs/${OUT}_\$TASK_ID.log
#\$ -q ${QUEUES}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
source ${CONDA}; conda activate ${ENV}
cd ${HERE}
python zak_sweep.py --index \$((SGE_TASK_ID - 1)) --grid ${GRID} --outdir ${OUT}
EOF

cat > script-plot-${OUT} <<EOF
#!/bin/bash
#\$ -clear
#\$ -S /bin/bash
#\$ -N zak_plot_${OUT}
#\$ -cwd
#\$ -l mem=2G
#\$ -j y
#\$ -o logs/plot_${OUT}.log
#\$ -q ${QUEUES}
#\$ -hold_jid zak_${OUT}
export OMP_NUM_THREADS=1 MPLBACKEND=Agg
source ${CONDA}; conda activate ${ENV}
cd ${HERE}
[ -d results ] && python migrate_old.py || true
python make_plots.py --results results_chiC        --grid grid_chiC.json        --plotdir plots_chiC        || true
python make_plots.py --results results_omegar      --grid grid_omegar.json      --plotdir plots_omegar      || true
python make_plots.py --results results_gammaomegar --grid grid_gammaomegar.json --plotdir plots_gammaomegar || true
EOF

echo "Array: ${NSWEEPS} tasks  grid=${GRID}  out=${OUT}  queues=${QUEUES}"
qsub script-array-${OUT}
qsub script-plot-${OUT}
echo "Watch: qstat -u \$USER ;  done when  ls ${OUT}/*/data.npz | wc -l  == ${NSWEEPS}"
