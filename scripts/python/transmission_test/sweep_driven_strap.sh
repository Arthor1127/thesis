#!/bin/bash
#$ -cwd
#$ -N drv_sweep
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=32G
#$ -pe neworte 16
#$ -t 1-200
#$ -tc 4
#$ -V

# Driven S21 sweep over BOTH axes:
#   - frequency: fixed 14-16 GHz band, same for every task
#   - strap position: one position per array task, centred on CENTER_MM
#     and spread +-AMPLITUDE_MM over N_STEPS values
#
# Positions are linspace(CENTER_MM - AMPLITUDE_MM, CENTER_MM + AMPLITUDE_MM,
# N_STEPS), so with N_STEPS odd the centre value is hit exactly.
#
# IMPORTANT: N_STEPS must match the "#$ -t 1-N" range above. SGE parses
# that directive before this script runs, so it cannot be derived from
# N_STEPS automatically - change both together.

CENTER_MM=3.9527        # 2L/3 - where both segments land on 3*f1 = 15 GHz
AMPLITUDE_MM=0.3        # sweep +- this far around the centre
N_STEPS=200             # MUST equal the -t range above

FREQ_START_GHZ=14.0
FREQ_END_GHZ=16.0

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg

PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/transmission_test"
cd "$SCRIPT_DIR"

# Position for THIS task. Validated against build_design.py's own
# EDGE_MARGIN_MM and total length rather than hardcoded here, so a bad
# CENTER/AMPLITUDE combination fails loudly at task start instead of
# raising ValueError deep inside build_design().
POS=$(python3 -c "
import sys
import numpy as np
from build_design import TOTAL_LENGTH_MM, EDGE_MARGIN_MM
center, amp, n = $CENTER_MM, $AMPLITUDE_MM, $N_STEPS
tid = $SGE_TASK_ID
if not (1 <= tid <= n):
    sys.exit(f'task id {tid} outside 1..{n} - does -t match N_STEPS?')
pos = center if n == 1 else float(np.linspace(center-amp, center+amp, n)[tid-1])
lo, hi = EDGE_MARGIN_MM, TOTAL_LENGTH_MM - EDGE_MARGIN_MM
if not (lo <= pos <= hi):
    sys.exit(f'position {pos:.5f} outside valid [{lo}, {hi:.5f}] - '
             f'reduce AMPLITUDE_MM or move CENTER_MM')
print(f'{pos:.5f}')
") || {
    echo "FATAL: could not compute strap position for task $SGE_TASK_ID" >&2
    exit 1
}

OUTDIR=$HOME/sweep_driven/pos_${POS}
mkdir -p $OUTDIR

echo "=== Driven task $SGE_TASK_ID/$N_STEPS: strap at ${POS}mm (from SHORT end) ==="
echo "=== Band ${FREQ_START_GHZ}-${FREQ_END_GHZ} GHz, centre=${CENTER_MM} amp=+-${AMPLITUDE_MM} ==="
echo "=== Running on host: $(hostname) ==="

python build_and_driven.py \
    --ground-pos $POS \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --freq-start-ghz $FREQ_START_GHZ \
    --freq-end-ghz $FREQ_END_GHZ \
    --n-points 2001 \
    --solver-order 2 \
    --adaptive-tol 1e-3 \
    --adaptive-max-samples 40 \
    --bg-min-size-um 8 --bg-max-size-um 100 \
    --resonator-bg-min-size-um 15 --resonator-bg-max-size-um 150 \
    --strap-taper-dist-max-um 40 \
    --solns-to-save 0

# Clean up heavy simulation files, keeping only the datasets
echo "Cleaning up heavy mesh and field files in $OUTDIR..."
find "$OUTDIR" -type f \( -name "*.msh" -o -name "*.pvu" -o -name "*.vtu" \) -delete
find "$OUTDIR" -type d -name "paraview" -exec rm -rf {} +