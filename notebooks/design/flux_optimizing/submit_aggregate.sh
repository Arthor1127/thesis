#!/bin/bash
#$ -cwd
#$ -N flux_aggregate
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=8G
#$ -pe neworte 2
#$ -V

# Runs after the sweep chain (via -hold_jid). Merges sweep/shard_*.csv ->
# sweep_results.csv, prints top configs + 1-D marginals, renders the 10 pairwise
# colormaps, and writes sweep/DONE.txt.

set -uo pipefail
conda deactivate 2>/dev/null || true
conda activate qcg-quantum-design
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PY=$HOME/.conda/envs/qcg-quantum-design/bin/python

cd "$HOME/circuit_design/flux_optimizing"
NEED=$($PY -c "import sys;sys.path.insert(0,'.');import sweep_flux;print(len(sweep_flux.all_points()))")
read DONE ERR < <($PY count_done.py sweep)

echo "=== aggregate $(date): done=$DONE/$NEED  errors=$ERR ===" | tee sweep/AGGREGATE.log
$PY aggregate_sweep.py --sweep-dir sweep 2>&1 | tee -a sweep/AGGREGATE.log
STATUS=${PIPESTATUS[0]}

{
  echo "flux sweep finished: $(date)"
  echo "grid points: $DONE / $NEED computed   (errors/missing: $((NEED - DONE)))"
  echo
  echo "  sweep/sweep_results.csv     - full table (last row per grid point)"
  echo "  sweep/sweep_colormaps.png   - 10 pairwise max-projected colormaps"
  echo "  sweep/AGGREGATE.log         - top-15 configs, bottom-5, 1-D marginal optima"
  echo
  echo "aggregate exit status: $STATUS"
} > sweep/DONE.txt
echo "=== aggregate done $(date) status=$STATUS ==="
