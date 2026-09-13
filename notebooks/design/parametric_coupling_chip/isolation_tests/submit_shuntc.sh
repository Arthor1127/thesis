#!/bin/bash
#$ -N shuntc
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=24G
#$ -pe neworte 16
#$ -o shuntc_task_$TASK_ID.log
#$ -t 1-4

# Validates what create_port_2_conds's impedance_C actually realizes -- see
# shuntc_validation_design.py for the full motivation (the two-port hanger
# isolation test's Qc came out wildly inconsistent vs its analytic target,
# and this isolates whether the port's capacitance realization is the cause,
# with zero resonator/mode-finding confounds). No resonance -> cheap,
# broadband, no PROM tuning needed.

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
mkdir -p out_shuntc

PARAMS_FILE="${PARAMS_FILE:-isolation_tests/params_shuntc.txt}"
LINE=$(sed -n "${SGE_TASK_ID}p" "$PARAMS_FILE")
read -r CC TAG PR <<< "$LINE"
export CC_FF="$CC" TAG PALACE_BIN
export PORT_R_OHM="${PR:-0}"
export FSTART_GHZ=1.0 FEND_GHZ=15.0 NPTS=141
export OUTDIR="out_shuntc"

echo "=== shuntc task $SGE_TASK_ID : Cc=${CC}fF -> TAG=$TAG host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u isolation_tests/run_shuntc_validation.py 2>&1
ST=$?
echo "=== shuntc task $SGE_TASK_ID ($TAG) exit=$ST $(date) ==="
exit $ST
