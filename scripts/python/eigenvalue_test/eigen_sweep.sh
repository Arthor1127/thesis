#!/bin/bash
#$ -cwd
#$ -N eigen_sweep
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem,cpuINgpu,be_cpuINgpu,be_copahue
#$ -l memoria_a_usar=12G
#$ -pe neworte 8
#$ -t 1-250
#$ -V

# Eigenvalue calculation sweep:
#   - strap position: one position per array task, centred on CENTER_MM
#     and spread +-AMPLITUDE_MM over N_STEPS values
#   - starting-freq/target-freq-ghz: computed PER TASK from the fitted
#     two-segment model, not fixed. Each task's shift-invert search is
#     centred right on its own predicted mode, instead of a single fixed
#     value that only matches the mode near the centre position. The
#     fitted envelope spans ~12.87-17.47 GHz across the full sweep range -
#     a fixed target anywhere near that would be badly off for positions
#     at the other end, exactly the mismatch that produced a header-only
#     eig.csv earlier in this project.
#
# Positions are linspace(CENTER_MM - AMPLITUDE_MM, CENTER_MM + AMPLITUDE_MM,
# N_STEPS), so with N_STEPS odd the centre value is hit exactly.
#
# IMPORTANT: N_STEPS must match the "#$ -t 1-N" range above. SGE parses
# that directive before this script runs, so it cannot be derived from
# N_STEPS automatically - change both together.

CENTER_MM=3.9527        # 2L/3 - where both segments land on 3*f1 = 15 GHz
AMPLITUDE_MM=0.3        # sweep +- this far around the centre
N_STEPS=250             # MUST equal the -t range above

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg

PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/eigenvalues_test"
cd "$SCRIPT_DIR"

# Position AND target frequency for THIS task, from the fitted two-segment
# model. v = 4*L*f_bare/3, back-calculated from the measured degenerate
# mode at x0=2L/3 (15.0324GHz) - see the anticrossing analysis this
# session. CORR_A/CORR_B (-0.6%/-1.4%) are the real residuals left after
# the ideal 1/(2x0) and 1/(4*(L-x0)) geometry factors are divided out -
# NOT free fit parameters, they correct for the strap's finite inductance
# (A) and open-end fringing capacitance (B).
read -r POS TARGET_FREQ_GHZ <<< $(python3 -c "
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

L = TOTAL_LENGTH_MM
f_bare = 15.0324
v = 4.0 * L * f_bare / 3.0          # GHz*mm
CORR_A, CORR_B = 0.994, 0.986

xA, xB = pos, L - pos
freq_a = CORR_A * v / (2.0 * xA)     # segment A: short-short
freq_b = CORR_B * v / (4.0 * xB)     # segment B: short-open
target = min(freq_a, freq_b)

print(f'{pos:.5f} {target:.5f}')
") || {
    echo "FATAL: could not compute strap position/target for task $SGE_TASK_ID" >&2
    exit 1
}

OUTDIR=$HOME/sweep_eigenmode/pos_${POS}
mkdir -p $OUTDIR

echo "=== Eigenmode task $SGE_TASK_ID/$N_STEPS: strap at ${POS}mm (from SHORT end) ==="
echo "=== Predicted lowest mode (fitted model): ${TARGET_FREQ_GHZ} GHz - using as both starting-freq and target-freq-ghz ==="
echo "=== Running on host: $(hostname) ==="

python build_and_eigenmode.py \
    --ground-pos $POS \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --starting-freq ${TARGET_FREQ_GHZ}e9 \
    --target-freq-ghz $TARGET_FREQ_GHZ \
    --number-of-freqs 3 \
    --amr-max-its 1 \
    --bg-min-size-um 8 --bg-max-size-um 100 \
    --resonator-bg-min-size-um 15 --resonator-bg-max-size-um 150 \
    --strap-taper-dist-max-um 40

# Clean up heavy simulation files, keeping only the datasets
echo "Cleaning up heavy mesh and field files in $OUTDIR..."
find "$OUTDIR" -type f \( -name "*.msh" -o -name "*.pvu" -o -name "*.vtu" \) -delete
find "$OUTDIR" -type d -name "paraview" -exec rm -rf {} +
