#!/bin/bash
#$ -N cpwcap_len
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
#$ -l memoria_a_usar=16G
#$ -l h_rss=12G
#$ -l h_rt=3600
#$ -pe neworte 8
#$ -tc 6
#$ -o cpwcap_len_task_$TASK_ID.log
#$ -t 1-6

# -o is a FLAT path on purpose (see submit_capacitance.sh's comment): SGE
# opens it before the job body's `mkdir -p` runs.

# Length-linearity check (capacitance side): same (width, gap) point at
# several trace lengths -- see gen_lengthcheck_params.py for the grid and
# check_length_linearity.py for the analysis. -t must match that script's
# printed row count (default 6 = 2 points x 3 lengths).

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
mkdir -p out_cap_lengthcheck

ulimit -v $((24 * 1024 * 1024)) || true

LINE=$(sed -n "${SGE_TASK_ID}p" params_lengthcheck.txt)
read -r W G LEN TAG <<< "$LINE"
export CPW_WIDTH_UM="$W" CPW_GAP_UM="$G" LENGTH_UM="$LEN" TAG PALACE_BIN
export OUTDIR="out_cap_lengthcheck"

( while true; do printf "%s %s\n" "$(date +%s)" \
    "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
    sleep 5; done ) > "out_cap_lengthcheck/ram_${SGE_TASK_ID}.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

echo "=== task $SGE_TASK_ID : $LINE -> TAG=$TAG  host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u run_capacitance.py 2>&1
ST=$?
kill $MON 2>/dev/null
echo "=== task $SGE_TASK_ID exit=$ST  peak RSS: $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' out_cap_lengthcheck/ram_${SGE_TASK_ID}.log 2>/dev/null) ($(date)) ==="
exit $ST
