#!/bin/bash
#$ -cwd
#$ -N lprobe_sweep_l600
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=96G
#$ -pe neworte 16
#$ -t 1-3
#$ -tc 1
#$ -V

# L_probe sweep at FIXED (manually adjusted) mesh - readout taper widened,
# squid-loop-specific fine mesh added, per the 2026-08-25 AMR findings -
# to check whether the extracted c parameter is stable across L_probe or
# shows loop-inductance-screening dependence (EPR_C_EXTRACTION.md Sec
# 6.1). No AMR here - single mesh+solve per task, comparable cost to one
# AMR iteration (~80G peak observed at similar mesh size), hence 96G.
# -tc 1 (strictly sequential, not parallel) since it's still daytime and
# other users may need the cluster - see 2026-08-25 chat.

set -uo pipefail

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_test"
cd "$SCRIPT_DIR"

LPROBE_VALUES=(10 15 20)
LPROBE=${LPROBE_VALUES[$((SGE_TASK_ID-1))]}
OUTDIR="$SCRIPT_DIR/lprobe_sweep_l1_600um/Lprobe_${LPROBE}nH"
mkdir -p "$OUTDIR"

echo "=== Task $SGE_TASK_ID: L_probe=${LPROBE}nH, host=$(hostname) ==="

timeout 3600 python -u build_and_eigenmode.py \
    --L-probe-nH "$LPROBE" \
    --outdir "$OUTDIR" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    --number-of-freqs 3 \
    --starting-freq 7.0e9 \
    2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}
echo "=== Task $SGE_TASK_ID done, exit=$STATUS ==="
