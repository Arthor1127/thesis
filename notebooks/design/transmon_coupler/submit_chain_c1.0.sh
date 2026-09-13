#!/bin/bash
#$ -N chainc10
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=110G
#$ -l h_rss=90G
#$ -l h_rt=21600
#$ -pe neworte 16
#$ -o chain_cap/task2_2.log
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

# The 1.0 um coupler-mesh point of the finer convergence check, split onto the
# copahue queue (16 slots) so it runs in parallel with the 1.5 um point on
# compute-4-18 instead of queueing behind it. Electrostatic solve scales fine
# to 16 ranks; RAM was ~14 GB at 1.5 M DoF so 90 GB h_rss is ample headroom.
export N_TRANSMONS=6
export COUPLER_RES_UM=1.0
export TRANSMON_RES_UM=12.0
export SOLVER_ORDER=2
export BG_MAX_UM=80
export MARGIN_X_MM=1.2 MARGIN_Y_MM=1.2
export TAG="chain6_c1.0_o2"

( while true; do printf "%s %s\n" "$(date +%s)" \
   "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
   sleep 10; done ) > "chain_cap/ram2_2.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

echo "=== $TAG  coupler_res=1.0um order 2 (copahue 16 slots)  $(date) ==="
python -u run_chain_cap.py 2>&1
echo "=== exit=$?  peak RSS(mon): $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' chain_cap/ram2_2.log) $(date) ==="
