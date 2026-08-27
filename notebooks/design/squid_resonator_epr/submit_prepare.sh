#!/bin/bash
#$ -cwd
#$ -N prep_unknowns
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=64G
#$ -pe neworte 16
#$ -V

# Goal: mesh + write the Palace config for the "full calculation" (finer
# mesh from run_coarse_test.py's edited parameters), then start Palace
# just long enough for it to print "number of global unknowns" (happens
# right after FE space assembly, well before the eigensolve iterations),
# then kill it. This tells us the problem size WITHOUT paying for the
# full (potentially very long) eigensolve.
#
# Memory: 64G, matching this account's own stageA_eigenmode.sh convention
# - NOT the ~450G requested in an earlier (overcorrected) attempt, which
# would have monopolized the cluster's only highmem node for a step that
# only OOM'd on a 14GB local workstation. Escalate only if this actually
# OOMs here too. Queue copahue,highmem (not highmem alone) so SGE can use
# any of the 6 copahue nodes rather than only compute-4-18.

set -uo pipefail

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg

PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"
export PALACE_BIN

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_test"
cd "$SCRIPT_DIR"

echo "=== Running on host: $(hostname) ==="
echo "=== NSLOTS: ${NSLOTS:-unset} ==="
echo "=== PALACE_BIN: $PALACE_BIN ==="

PREP_OUTDIR="$SCRIPT_DIR/prep_run"
mkdir -p "$PREP_OUTDIR"

echo "=== [A] Building design + generating mesh + writing Palace config (no solve) ==="
python -u prepare_full_run.py "$PREP_OUTDIR" 2>&1 | tee "$PREP_OUTDIR/prep.log"
PREP_STATUS=${PIPESTATUS[0]}
if [ "$PREP_STATUS" -ne 0 ]; then
    echo "!!! prepare_full_run.py failed with exit code $PREP_STATUS - stopping here."
    exit 1
fi

SIM_DIR="$PREP_OUTDIR/full_run1_prep"
CONFIG_JSON="$SIM_DIR/full_run1_prep.json"
if [ ! -f "$CONFIG_JSON" ]; then
    echo "!!! Expected config not found at $CONFIG_JSON - listing $SIM_DIR:"
    ls -la "$SIM_DIR" 2>&1
    exit 1
fi

echo "=== [B] Mesh file size: ==="
ls -la "$SIM_DIR"/*.msh

echo "=== [C] Starting Palace just to get 'number of global unknowns', then killing it ==="
cd "$SIM_DIR"
PALACE_LOG="$SIM_DIR/palace_unknowns_probe.log"
"$PALACE_BIN" -np "${NSLOTS:-16}" full_run1_prep.json > "$PALACE_LOG" 2>&1 &
PALACE_PID=$!

FOUND=0
for i in $(seq 1 360); do   # up to 30 min (5s * 360), mesh assembly should be much faster
    if grep -q "number of global unknowns" "$PALACE_LOG" 2>/dev/null; then
        sleep 3   # let it print the "Level 0 (p = N): X unknowns" follow-up lines too
        FOUND=1
        break
    fi
    if ! kill -0 "$PALACE_PID" 2>/dev/null; then
        echo "!!! Palace process exited before printing unknown count."
        break
    fi
    sleep 5
done

echo "=== Killing Palace/mpirun (process group) ==="
kill -TERM -"$PALACE_PID" 2>/dev/null
sleep 2
kill -KILL -"$PALACE_PID" 2>/dev/null
pkill -9 -f "full_run1_prep.json" 2>/dev/null || true

echo "=== Palace probe log (relevant lines) ==="
grep -iE "unknowns|error|Assembling|memory|Mesh|MPI|Backend" "$PALACE_LOG" | head -60

if [ "$FOUND" -eq 1 ]; then
    echo "=== SUCCESS: unknown count captured above ==="
else
    echo "=== WARNING: did not find 'number of global unknowns' within timeout - check $PALACE_LOG in full ==="
fi
