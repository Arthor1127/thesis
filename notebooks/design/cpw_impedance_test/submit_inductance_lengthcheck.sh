#!/bin/bash
#$ -N cpwind_len
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
# (was copahue -- saturated by the flux_sweep campaign; highmem has slack and
# already ran the capacitance jobs fine)
#$ -l memoria_a_usar=16G
#$ -pe neworte 8
#$ -tc 6
#$ -o cpwind_len_task_$TASK_ID.log
#$ -t 1-6

# -o is a FLAT path on purpose (see submit_capacitance.sh's comment): SGE
# opens it before the job body's `mkdir -p` runs.

# Length-linearity check (inductance side): same (width, gap) point at
# several trace lengths -- see gen_lengthcheck_params.py for the grid and
# check_length_linearity.py for the analysis. -t must match that script's
# printed row count (default 6 = 2 points x 3 lengths). Note the flux
# integration area is already margin-trimmed per length by
# build_design.integration_area_gap_strips, so integration_length_um (not
# length_um) is the fair x-axis for the L_total-vs-length fit -- see
# check_length_linearity.py.

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
mkdir -p out_ind_lengthcheck

LINE=$(sed -n "${SGE_TASK_ID}p" params_lengthcheck.txt)
read -r W G LEN TAG <<< "$LINE"

( while true; do printf '%s ' "$(date +%H:%M:%S)"; free -g | awk 'NR==2{printf "used=%sG avail=%sG\n",$3,$7}'; sleep 30; done ) > "out_ind_lengthcheck/ram_${SGE_TASK_ID}.log" 2>&1 &
RAMPID=$!; trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== task $SGE_TASK_ID : $LINE -> TAG=$TAG  host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 7200 python -u run_inductance.py \
    --outdir out_ind_lengthcheck \
    --palace-bin "$PALACE_BIN" \
    --cpw-width-um "$W" \
    --cpw-gap-um "$G" \
    --length-um "$LEN" \
    --tag "$TAG" \
    --num-cpus "$NSLOTS" \
    2>&1 | tee "out_ind_lengthcheck/run_${SGE_TASK_ID}.log"
STATUS=${PIPESTATUS[0]}
kill $RAMPID 2>/dev/null
echo "=== task $SGE_TASK_ID done exit=$STATUS peakRAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m"G"}' out_ind_lengthcheck/ram_${SGE_TASK_ID}.log) $(date) ==="
exit $STATUS
