#!/bin/bash
#$ -N chainamr
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
#$ -l memoria_a_usar=180G
#$ -l h_rss=150G
#$ -l h_rt=36000
#$ -pe neworte 32
#$ -o chain_cap/amr.log
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

export N_TRANSMONS=6
export COUPLER_RES_UM=2.5
export TRANSMON_RES_UM=12.0
export SOLVER_ORDER=2
export BG_MAX_UM=80
export MARGIN_X_MM=1.2 MARGIN_Y_MM=1.2
export AMR_ITERS=6
export AMR_TOL=5e-3
export AMR_DORFLER=0.7
export AMR_MAX_DOFS=12000000       # RAM guard: ~12M order-2 DoFs ~ <120 GB
export TAG="chain6_amr_o2"

( while true; do printf "%s %s\n" "$(date +%s)" \
   "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
   sleep 15; done ) > "chain_cap/ram_amr.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

echo "=== $TAG  AMR 6 iters tol 5e-3  $(date) ==="
python -u run_chain_cap_amr.py 2>&1
echo "=== exit=$?  peak RSS(mon): $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' chain_cap/ram_amr.log) $(date) ==="
