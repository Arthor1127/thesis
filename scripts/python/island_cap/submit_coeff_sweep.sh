#!/bin/bash
#$ -N islcoef
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
#$ -l memoria_a_usar=60G
#$ -l h_rss=48G
#$ -l h_rt=2400
#$ -pe neworte 32
#$ -tc 1
#$ -o coeff/task_$TASK_ID.log
set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
export PALACE_BIN="$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace"
export OUTDIR="coeff"
cd ~/circuit_design/island_cap
mkdir -p coeff/results

# Refine the coefficients of  C/(eps_eff*a) = A + B*ln(a/g)  from ORDER-2,
# mesh-converged runs (validated by the a=136/g=20 cross-check: 97.0 fF vs the
# formula's 100 fF, stable 2 vs 3 um mesh). Same mesh recipe: order 2, thin
# rim ring, res_min = max(g/6, 1.2) um  (~6 elements across the gap).
LINE=$(sed -n "${SGE_TASK_ID}p" params_coeff.txt)
read -r A G <<< "$LINE"
export RADIUS_UM="$A" GAP_UM="$G" THICK_NM=100 SOLVER_ORDER=2
export RES_MIN_UM=$(python3 -c "print(max($G/6.0, 1.2))")
export TAG="coeff_a${A}_g${G}_o2"

( while true; do printf "%s %s\n" "$(date +%s)" \
   "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
   sleep 5; done ) > "coeff/ram_${SGE_TASK_ID}.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

echo "=== $TAG  res_min=${RES_MIN_UM}um order 2  $(date) ==="
python -u run_cap_island.py 2>&1
echo "=== task $SGE_TASK_ID exit=$?  peak RSS: $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' coeff/ram_${SGE_TASK_ID}.log) $(date) ==="
