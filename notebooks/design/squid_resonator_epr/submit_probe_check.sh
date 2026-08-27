#!/bin/bash
#$ -cwd
#$ -N probe_check
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=64G
#$ -pe neworte 16
#$ -V
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
OUTDIR="$SCRIPT_DIR/probe_check"
mkdir -p "$OUTDIR"

# Run in background so we can kill it right after the unknown count prints,
# before the expensive multigrid/solve phase - avoids repeating the near-OOM.
timeout 400 python -u build_and_eigenmode.py \
    --L-probe-nH 15 --outdir "$OUTDIR" --palace-bin "$PALACE_BIN" --num-cpus "$NSLOTS" \
    --number-of-freqs 3 --starting-freq 7.0e9 \
    > "$OUTDIR/run.log" 2>&1 &
PYPID=$!

FOUND=0
for i in $(seq 1 78); do
    if grep -q "number of global unknowns" "$OUTDIR/run.log" 2>/dev/null; then
        sleep 2
        FOUND=1
        break
    fi
    if ! kill -0 "$PYPID" 2>/dev/null; then
        echo "!!! process exited before printing unknown count"
        break
    fi
    sleep 5
done

echo "=== Result ==="
grep -A1 "number of global unknowns" "$OUTDIR/run.log" | tail -5

echo "=== Killing ==="
kill -TERM -"$PYPID" 2>/dev/null
sleep 2
kill -KILL -"$PYPID" 2>/dev/null
pkill -9 -f "palace-x86_64.bin" 2>/dev/null || true
pkill -9 -f "build_and_eigenmode.py" 2>/dev/null || true

if [ "$FOUND" -eq 0 ]; then
    echo "=== WARNING: never found unknown count ==="
    tail -n 30 "$OUTDIR/run.log"
fi
