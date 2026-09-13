#!/bin/bash
#$ -N chaincap2
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
#$ -l memoria_a_usar=110G
#$ -l h_rss=90G
#$ -l h_rt=21600
#$ -pe neworte 32
#$ -tc 1
#$ -o chain_cap/task2_$TASK_ID.log
set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
export PALACE_BIN="$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace"
export OUTDIR="chain_cap"
cd ~/circuit_design/transmon_coupler
mkdir -p chain_cap/results

# Finer mesh-convergence: coupler fine-mesh size 1.5 um then 1.0 um (the 3.0/2.0
# pair was NOT converged for the interdigital-coupler capacitances).
CPL_RES=$(sed -n "${SGE_TASK_ID}p" <<< $'1.5\n1.0')
export N_TRANSMONS=6
export COUPLER_RES_UM="$CPL_RES"
export TRANSMON_RES_UM=12.0
export SOLVER_ORDER=2
export BG_MAX_UM=80
export MARGIN_X_MM=1.2 MARGIN_Y_MM=1.2
export TAG="chain6_c${CPL_RES}_o2"

( while true; do printf "%s %s\n" "$(date +%s)" \
   "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
   sleep 10; done ) > "chain_cap/ram2_${SGE_TASK_ID}.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

echo "=== $TAG  coupler_res=${CPL_RES}um order 2  $(date) ==="
python -u run_chain_cap.py 2>&1
echo "=== task $SGE_TASK_ID exit=$?  peak RSS(mon): $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' chain_cap/ram2_${SGE_TASK_ID}.log) $(date) ==="
