#!/bin/bash
#$ -N cpwcap
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
#$ -l memoria_a_usar=16G
#$ -l h_rss=12G
#$ -l h_rt=3600
#$ -pe neworte 8
#$ -tc 6
#$ -o cpwcap_task_$TASK_ID.log
#$ -t 1-30

# -o is a FLAT path (not out_cap/task_...) on purpose: SGE opens it before the
# job body runs, so a subdirectory that only this script's own `mkdir -p`
# creates does not exist yet -> job goes to Eqw with "No such file or
# directory". Learned the hard way on the first submission of this script.

# CPW capacitance-per-length sweep over (cpw_width_um, cpw_gap_um) -- see
# gen_params.py for the grid and run_capacitance.py for the sim itself.
# Resource sizing is a first guess (~8 slots, well under the island_cap
# precedent's order-1 "peaks well under 10 GB" for a similarly simple
# 2-conductor electrostatic problem) -- pilot task 1 alone before
# launching the full -t 1-30 array, and tighten memoria_a_usar/h_rt from
# what it actually reports. h_rt was doubled from the pre-AMR version since
# AMR (on by default in run_capacitance.py) means up to AMR_ITERS mesh+solve
# passes per task, not one.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/cpw_impedance_test"
cd "$SCRIPT_DIR"
mkdir -p out_cap

# Loose virtual-memory backstop, same spirit as island_cap/submit_island_array.sh:
# only catches a genuine runaway, normal tasks here should stay far under it.
ulimit -v $((24 * 1024 * 1024)) || true

LINE=$(sed -n "${SGE_TASK_ID}p" params.txt)
read -r W G TAG <<< "$LINE"
export CPW_WIDTH_UM="$W" CPW_GAP_UM="$G" TAG PALACE_BIN
export OUTDIR="out_cap"

( while true; do printf "%s %s\n" "$(date +%s)" \
    "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
    sleep 5; done ) > "out_cap/ram_${SGE_TASK_ID}.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

echo "=== task $SGE_TASK_ID : $LINE -> TAG=$TAG  host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u run_capacitance.py 2>&1
ST=$?
kill $MON 2>/dev/null
echo "=== task $SGE_TASK_ID exit=$ST  peak RSS: $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' out_cap/ram_${SGE_TASK_ID}.log 2>/dev/null) ($(date)) ==="
exit $ST
