#!/bin/bash
#$ -cwd
#$ -N amr_run2
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=150G
#$ -pe neworte 16
#$ -V

# AMR run #2: readout taper widened (40um) / qubit taper narrowed (8um)
# per the 2026-08-25 AMR-1 Indicator-field findings. Capped at MaxIts=2
# (not 3) and deliberately kept off the highmem queue - it's daytime and
# other users may need the cluster, so this stays on copahue (multiple
# nodes with 141-220G) rather than monopolizing the single highmem node.
# 150G requested based on AMR-1's measured scaling (iter1 11.8G -> iter2
# 38.2G -> iter3 160G for a SMALLER base mesh; this run's base mesh is
# already 2.47M elements vs. iter1's 434k, so iteration 2 here could
# plausibly need well over 100G).
#
# Being watched live this time (not an unattended overnight run), so if
# it looks like it's headed toward exceeding the request or straining
# the node, it'll be killed manually rather than left to run.

set -uo pipefail

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
export PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"
export AMR_ITERS=2

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_test"
cd "$SCRIPT_DIR"

echo "=== Running on host: $(hostname) ==="
echo "=== NSLOTS: ${NSLOTS:-unset} ==="
echo "=== Start time: $(date) ==="

OUTDIR="$SCRIPT_DIR/amr_run2"
mkdir -p "$OUTDIR"

# 2-hour cap (shorter than the overnight run's 3h - this is being watched
# live, so if it's not done by then something is wrong and should be
# investigated rather than left running).
timeout 7200 python -u run_amr.py "$OUTDIR" 2>&1 | tee "$OUTDIR/amr.log"
STATUS=${PIPESTATUS[0]}

echo "=== End time: $(date) ==="
if [ "$STATUS" -eq 124 ]; then
    echo "!!! TIMED OUT after 2 hours"
elif [ "$STATUS" -ne 0 ]; then
    echo "!!! FAILED with exit code $STATUS"
else
    echo "=== SUCCESS ==="
fi
