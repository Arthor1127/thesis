#!/bin/bash
#$ -N qc_sweep
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue,highmem
#$ -l memoria_a_usar=32G
#$ -pe neworte 16
#$ -tc 4
#$ -o qc_task_$TASK_ID.log
#$ -t 1-35

# Qc sweep: driven S21 sim of the hanger resonator over (coupling_space_um,
# coupling_length_um) -- see gen_qc_sweep_params.py for the grid and
# run_qc_sim.py for the sim itself (fit_transmission gives fr, Qr, Qc, Qi).
#
# -o is a FLAT path on purpose: SGE opens it before the job body's `mkdir -p`
# runs (see notebooks/design/cpw_impedance_test/submit_capacitance.sh's
# comment for the Eqw failure this avoids).
#
# Resource sizing matches scripts/python/transmission_test/sweep_driven_strap.sh
# (32G / 16 slots for a driven sim) -- a first guess for THIS geometry, though;
# pilot task 1 alone before the full -t 1-35 array. -q lists both copahue and
# highmem so SGE can use whichever has free slots (copahue can be fully
# saturated by other campaigns -- see the cpw_impedance_test session notes).

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_chip"
cd "$SCRIPT_DIR"
mkdir -p out_qc

LINE=$(sed -n "${SGE_TASK_ID}p" driven_sweep/params_qc.txt)
read -r CS CL TAG <<< "$LINE"
export COUPLING_SPACE_UM="$CS" COUPLING_LENGTH_UM="$CL" TAG PALACE_BIN
export OUTDIR="out_qc"

( while true; do printf "%s %s\n" "$(date +%s)" \
    "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
    sleep 5; done ) > "out_qc/ram_${SGE_TASK_ID}.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

echo "=== task $SGE_TASK_ID : $LINE -> TAG=$TAG  host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u driven_sweep/run_qc_sim.py 2>&1
ST=$?
kill $MON 2>/dev/null
echo "=== task $SGE_TASK_ID exit=$ST  peak RSS: $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' out_qc/ram_${SGE_TASK_ID}.log 2>/dev/null) ($(date)) ==="
exit $ST
