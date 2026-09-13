#!/bin/bash
#$ -N cpwind
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
# (was copahue -- saturated by the flux_sweep campaign; highmem has slack and
# already ran the capacitance jobs fine)
#$ -l memoria_a_usar=16G
#$ -pe neworte 8
#$ -tc 6
#$ -o cpwind_task_$TASK_ID.log
#$ -t 1-30

# -o is a FLAT path on purpose: SGE opens it before the job body's `mkdir -p`
# runs, so a not-yet-existing subdirectory (out_ind/task_...) puts the job in
# Eqw with "No such file or directory". Learned the hard way on the first
# submission of this script.

# CPW self-inductance-per-length sweep over (cpw_width_um, cpw_gap_um) --
# same grid as submit_capacitance.sh (gen_params.py / params.txt), so a given
# SGE_TASK_ID is the same (width, gap) point in both sweeps. See
# run_inductance.py for the method (flux linkage through the CPW gap, current
# driven at the open end, returned via the short).
#
# Resource sizing is a first guess -- a plain straight-line loop (no
# topological handles, unlike flux_optimizing's SQUID array) should mesh/solve
# markedly lighter than that precedent's 32-slot/32G runs; pilot task 1 alone
# before the full -t 1-30 array and tighten from what it actually reports.
# The `timeout` below was doubled from the pre-AMR version since AMR (on by
# default in run_inductance.py) means multiple mesh+solve passes per task.

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
mkdir -p out_ind

LINE=$(sed -n "${SGE_TASK_ID}p" params.txt)
read -r W G TAG <<< "$LINE"

( while true; do printf '%s ' "$(date +%H:%M:%S)"; free -g | awk 'NR==2{printf "used=%sG avail=%sG\n",$3,$7}'; sleep 30; done ) > "out_ind/ram_${SGE_TASK_ID}.log" 2>&1 &
RAMPID=$!; trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== task $SGE_TASK_ID : $LINE -> TAG=$TAG  host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 7200 python -u run_inductance.py \
    --outdir out_ind \
    --palace-bin "$PALACE_BIN" \
    --cpw-width-um "$W" \
    --cpw-gap-um "$G" \
    --tag "$TAG" \
    --num-cpus "$NSLOTS" \
    2>&1 | tee "out_ind/run_${SGE_TASK_ID}.log"
STATUS=${PIPESTATUS[0]}
kill $RAMPID 2>/dev/null
echo "=== task $SGE_TASK_ID done exit=$STATUS peakRAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m"G"}' out_ind/ram_${SGE_TASK_ID}.log) $(date) ==="
exit $STATUS
