#!/bin/bash
#$ -cwd
#$ -N sqarr_ljscan
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=48G
#$ -pe neworte 32
#$ -t 1-3
#$ -tc 3
#$ -V

# Coarse, no-AMR scan over junction inductance L_J to find a value where the
# whole SQUID-array phonon band sits BELOW the readout resonator (~8 GHz on
# the refined mesh), so the resonator mode doesn't hybridise with a phonon.
# Phonon band top ~ 1/sqrt(L_J), so larger L_J pushes it down. Solves a
# window from 4.5 GHz (catches the top of the phonon band + the resonator).

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"
cd "$SCRIPT_DIR"

LJ_VALUES=(250 400 600)
LJ=${LJ_VALUES[$((SGE_TASK_ID-1))]}
OUTDIR="$SCRIPT_DIR/ljscan/Lj_${LJ}nH"
mkdir -p "$OUTDIR"

echo "=== Task $SGE_TASK_ID: L_J=${LJ}nH, host=$(hostname), $(date) ==="

timeout 2400 python -u build_and_eigenmode.py \
    --L-j-nH "$LJ" \
    --outdir "$OUTDIR" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    --number-of-freqs 8 \
    --starting-freq 4.5e9 \
    --bg-min-size-um 8 \
    --resonator-bg-min-size-um 8 \
    --array-bg-min-size-um 6 \
    --fillet-resolution 8 \
    2>&1 | tee "$OUTDIR/run.log"
echo "=== Task $SGE_TASK_ID done, exit=${PIPESTATUS[0]}, $(date) ==="
