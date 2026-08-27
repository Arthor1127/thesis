#!/bin/bash
#$ -cwd
#$ -N combined_6um
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=64G
#$ -pe neworte 16
#$ -V

# All four fine-mesh regions (readout, pad, squid_arm, junction) at the
# 6um floor identified by the 2026-08-25 per-region sweep (each region
# individually passed at 6um in 80-96s / ~5.5-6.2M elements, all timed
# out at 5um). This combined run will be larger than any single-region
# test since none of them were run together before - hence the generous
# (not 2-min) timeout on meshing, and 64G mem (same right-sized figure
# used successfully before, not escalated without evidence).
#
# If meshing succeeds, briefly starts Palace just long enough to print
# "number of global unknowns", then kills it - same probe technique as
# the 2026-08-24 single-run investigation, so we get problem size before
# committing to the full (possibly long) eigensolve.

set -uo pipefail

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
export PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_test"
cd "$SCRIPT_DIR"

echo "=== Running on host: $(hostname) ==="
echo "=== NSLOTS: ${NSLOTS:-unset} ==="

OUTDIR="$SCRIPT_DIR/combined_6um_run"
mkdir -p "$OUTDIR"

echo "=== [A] Mesh generation, all 4 regions at 6um (timeout 600s) ==="
timeout 600 python -u mesh_sweep_step.py --all-min-size-um 6 --outdir "$OUTDIR" 2>&1 | tee "$OUTDIR/mesh.log"
MESH_STATUS=${PIPESTATUS[0]}

if [ "$MESH_STATUS" -eq 124 ]; then
    echo "!!! Mesh generation TIMED OUT after 600s - stopping here."
    exit 1
elif [ "$MESH_STATUS" -ne 0 ]; then
    echo "!!! Mesh generation FAILED (exit $MESH_STATUS) - stopping here."
    exit 1
fi

SIM_DIR="$OUTDIR/sweep_step"
CONFIG_JSON="$SIM_DIR/sweep_step.json"
if [ ! -f "$CONFIG_JSON" ]; then
    echo "!!! Expected config not found at $CONFIG_JSON"
    ls -la "$SIM_DIR" 2>&1
    exit 1
fi
echo "=== Mesh file: ==="
ls -la "$SIM_DIR"/*.msh

echo "=== [B] Starting Palace, letting it reach the iterative eigensolve (not just the unknown count), sampling live RSS, then killing it ==="
cd "$SIM_DIR"
PALACE_LOG="$SIM_DIR/palace_unknowns_probe.log"
MEM_TRACE="$SIM_DIR/rss_trace.log"
"$PALACE_BIN" -np "${NSLOTS:-16}" sweep_step.json > "$PALACE_LOG" 2>&1 &
PALACE_PID=$!

sample_rss_gb() {
    ps -eo rss,cmd | grep "palace-x86_64.bin" | grep -v grep | awk '{s+=$1} END {printf "%.2f", s/1024/1024}'
}

FOUND_UNKNOWNS=0
FOUND_SOLVE=0
SOLVE_START_ITER=0
SLEEP_S=5
WAIT_FOR_SOLVE_ITERS=96   # 96*5s = 480s cap waiting for the iterative solve to start
SAMPLE_ITERS=18           # 18*5s = 90s of RSS sampling once it has started
# Let it run past just the unknown-count line, into the (memory-heavy)
# multigrid hierarchy assembly and the start of the iterative eigensolve,
# sampling live RSS along the way, so we get a real measured peak instead
# of extrapolating from the small reference run - see 2026-08-25 chat.
for i in $(seq 1 $((WAIT_FOR_SOLVE_ITERS + SAMPLE_ITERS))); do
    if [ "$FOUND_UNKNOWNS" -eq 0 ] && grep -q "number of global unknowns" "$PALACE_LOG" 2>/dev/null; then
        FOUND_UNKNOWNS=1
        echo "[$(date +%T)] unknowns line seen, rss=$(sample_rss_gb)G" | tee -a "$MEM_TRACE"
    fi
    if [ "$FOUND_SOLVE" -eq 0 ] && grep -qE "Residual norms for FGMRES solve|KSP residual norm" "$PALACE_LOG" 2>/dev/null; then
        FOUND_SOLVE=1
        SOLVE_START_ITER=$i
        echo "[$(date +%T)] iterative solve started, rss=$(sample_rss_gb)G" | tee -a "$MEM_TRACE"
    fi
    if [ "$FOUND_SOLVE" -eq 1 ]; then
        elapsed_since_solve=$(( (i - SOLVE_START_ITER) * SLEEP_S ))
        echo "[$(date +%T)] +${elapsed_since_solve}s into solve, rss=$(sample_rss_gb)G" | tee -a "$MEM_TRACE"
        if [ $((i - SOLVE_START_ITER)) -ge "$SAMPLE_ITERS" ]; then
            break
        fi
    elif [ "$i" -ge "$WAIT_FOR_SOLVE_ITERS" ]; then
        echo "!!! Solve never started within $((WAIT_FOR_SOLVE_ITERS * SLEEP_S))s of Palace launch."
        break
    fi
    if ! kill -0 "$PALACE_PID" 2>/dev/null; then
        echo "!!! Palace process exited on its own (check log for error/OOM-kill)."
        break
    fi
    sleep "$SLEEP_S"
done

echo "=== Peak RSS observed while sampling: $(sort -t= -k2 -n "$MEM_TRACE" 2>/dev/null | tail -1) ==="
echo "=== Killing Palace/mpirun ==="
kill -TERM -"$PALACE_PID" 2>/dev/null
sleep 2
kill -KILL -"$PALACE_PID" 2>/dev/null
pkill -9 -f "sweep_step.json" 2>/dev/null || true
pkill -9 -f "palace-x86_64.bin" 2>/dev/null || true

echo "=== Palace probe log (relevant lines) ==="
grep -iE "unknowns|error|Assembling|memory|Mesh|MPI|Backend" "$PALACE_LOG" | head -60

echo "=== Full RSS trace ==="
cat "$MEM_TRACE" 2>/dev/null

if [ "$FOUND_UNKNOWNS" -eq 1 ] && [ "$FOUND_SOLVE" -eq 1 ]; then
    echo "=== SUCCESS: reached the iterative solve, RSS trace captured above ==="
elif [ "$FOUND_UNKNOWNS" -eq 1 ]; then
    echo "=== PARTIAL: got unknown count but never saw the iterative solve start ==="
else
    echo "=== WARNING: did not find 'number of global unknowns' within timeout ==="
fi
