#!/bin/bash
#$ -cwd
#$ -N jjarr_eig
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=180G
#$ -pe neworte 32
#$ -t 1-4
#$ -tc 2
#$ -V

# Normal modes (Palace eigenmode solve) of the bare Josephson-junction array
# (JosephsonArray component), sweeping the chain length. Every junction gap
# gets a lumped L_j = 1.2 nH / C_j = 2 fF port; the chain sits alone in a 50um
# ground-plane pocket. No resonator, no EPR / dielectric-loss analysis - just
# the eigenfrequencies (+ per-junction participation of each mode).
#
# task:        1     2     3     4
#  n_junctions:10    23    40    100
#
# N=100 is 100 lumped ports on a ~0.82mm chain; if it OOMs at 180G, re-run
# task 4 alone with -q highmem -l memoria_a_usar=460G.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

N_ARR=(10 23 40 100)
i=$((SGE_TASK_ID-1)); NJJ=${N_ARR[$i]}

SCRIPT_DIR="$HOME/circuit_design/josephson_array_resonator"; cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/eig/n${NJJ}"; mkdir -p "$OUTDIR"

( while true; do printf '%s ' "$(date +%H:%M:%S)"; free -g|awk 'NR==2{printf "used=%sG avail=%sG\n",$3,$7}'; sleep 30; done ) > "$OUTDIR/ram.log" 2>&1 &
RAMPID=$!; trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== jjarr_eig task $SGE_TASK_ID: n_junctions=${NJJ}, host=$(hostname), $(date) ==="
timeout 21600 python -u eigenmodes.py \
    --n-junctions "$NJJ" \
    --L-j-nH 1.2 --Cj-j-fF 2.0 \
    --ground-gap-um 50 \
    --number-of-freqs 12 --starting-freq 0.5e9 \
    --array-bg-min-size-um 2.0 --array-taper-dist-min-um 10 \
    --bg-min-size-um 8 --fillet-resolution 10 \
    --outdir "$OUTDIR" --palace-bin "$PALACE_BIN" --num-cpus "$NSLOTS" \
    2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}; kill $RAMPID 2>/dev/null
echo "=== task $SGE_TASK_ID done exit=$STATUS peakRAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m}' $OUTDIR/ram.log)G $(date) ==="
exit $STATUS
