#!/bin/bash
#$ -N islcap
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
#$ -l memoria_a_usar=20G
#$ -l h_rss=16G
#$ -l h_rt=1800
#$ -pe neworte 8
#$ -tc 3
#$ -o out/task_$TASK_ID.log
set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
export PALACE_BIN="$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace"
export OUTDIR="out"
cd ~/circuit_design/island_cap
mkdir -p out/results

# Loose virtual-memory backstop: only catches a genuine runaway (an order-2 run
# previously blew to 68 GB / 27M unknowns). Normal order-1 tasks here peak well
# under 10 GB. Runs on the 503 GB highmem node with -tc 4, so 4 * this is still
# a small fraction of the node -- the node cannot be hurt.
ulimit -v $((45 * 1024 * 1024)) || true

LINE=$(sed -n "${SGE_TASK_ID}p" params.txt)
read -r R G T SWEEP <<< "$LINE"
export RADIUS_UM="$R" GAP_UM="$G" THICK_NM="$T"
export TAG="task${SGE_TASK_ID}_r${R%.*}_g${G%.*}_t${T%.*}_${SWEEP}"

# lightweight RAM monitor for this task
( while true; do printf "%s %s\n" "$(date +%s)" \
    "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
    sleep 5; done ) > "out/ram_${SGE_TASK_ID}.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

echo "=== task $SGE_TASK_ID : $LINE  -> TAG=$TAG ($(date)) ==="
python -u run_cap_island.py 2>&1
ST=$?
kill $MON 2>/dev/null
echo "=== task $SGE_TASK_ID exit=$ST  peak RSS: $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' out/ram_${SGE_TASK_ID}.log 2>/dev/null) ($(date)) ==="
