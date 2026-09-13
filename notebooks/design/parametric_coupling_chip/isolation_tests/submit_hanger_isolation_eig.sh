#!/bin/bash
#$ -N hanger_eig
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=24G
#$ -pe neworte 16
#$ -o hanger_eig_task_$TASK_ID.log
#$ -t 1-3

# Two-port hanger-topology Qc isolation test, eigenmode pass -- see
# hanger_isolation_design.py and run_hanger_isolation_eigenmode.py for the
# full motivation. Run FIRST (before the driven-S21 companion) to pin down
# each point's real fr, needed to center a tight, cheap driven window.

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
mkdir -p out_hanger_eig

LINE=$(sed -n "${SGE_TASK_ID}p" isolation_tests/params_hanger_isolation_eig.txt)
read -r QC TAG <<< "$LINE"
export QC_TARGET="$QC" TAG PALACE_BIN
export STARTING_FREQ_GHZ=7.3 NUM_FREQS=2
export OUTDIR="out_hanger_eig"

echo "=== hanger_eig task $SGE_TASK_ID : Qc_target=$QC -> TAG=$TAG host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u isolation_tests/run_hanger_isolation_eigenmode.py 2>&1
ST=$?
echo "=== hanger_eig task $SGE_TASK_ID ($TAG) exit=$ST $(date) ==="
exit $ST
