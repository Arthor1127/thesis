#!/bin/bash
# ===========================================================================
#  Submit the stanene gap-map sweep array + a dependent plotting job.
#  Run ON the SGE cluster frontend, from  ~/unitary/cris/
#
#  Usage:  bash submit_gaps.sh ["queue1,queue2,..."] [tc]
#
#  1 slot/task, one task per B row (mesh tasks), checkpointed (resumable).
#  Verbatim notebook kernel; single-threaded BLAS (small dense linalg).
#  The `cris-gaps` conda env already exists on this cluster (pythtb stack).
# ===========================================================================
set -e
QUEUES=${1:-"cpuINgpu,copahue,be_cpuINgpu,be_copahue"}     # no highmem (other agent's)
TC=${2:-96}
ENV=cris-gaps
CONDA=/share/apps/miniconda3/etc/profile.d/conda.sh
OUT=results_stanene
HERE=$(pwd)

source ${CONDA}; conda activate ${ENV}
python build_grid.py
MESH=$(python -c "import json;print(json.load(open('grid.json'))['mesh'])")
mkdir -p logs ${OUT}

cat > script-array <<EOF
#!/bin/bash
#\$ -clear
#\$ -S /bin/bash
#\$ -N cris_stan
#\$ -cwd
#\$ -t 1-${MESH}
#\$ -tc ${TC}
#\$ -l mem=4G
#\$ -j y
#\$ -o logs/stan_row_\$TASK_ID.log
#\$ -q ${QUEUES}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
source ${CONDA}; conda activate ${ENV}
cd ${HERE}
python gaps_sweep.py --index \$((SGE_TASK_ID - 1)) --grid grid.json --outdir ${OUT}
EOF

cat > script-plot <<EOF
#!/bin/bash
#\$ -clear
#\$ -S /bin/bash
#\$ -N cris_stan_plot
#\$ -cwd
#\$ -l mem=2G
#\$ -j y
#\$ -o logs/stan_plot.log
#\$ -q ${QUEUES}
#\$ -hold_jid cris_stan
export OMP_NUM_THREADS=1 MPLBACKEND=Agg
source ${CONDA}; conda activate ${ENV}
cd ${HERE}
python make_plots.py --results ${OUT} --grid grid.json --outdir ${OUT}
echo "PLOT DONE"
EOF

echo "Array: ${MESH} row tasks  queues=${QUEUES}  tc=${TC}  out=${OUT}"
qsub script-array
qsub script-plot
echo "Watch: qstat -u \$USER ;  done when  ls ${OUT}/row_*.npz | wc -l  == ${MESH}"
