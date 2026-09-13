#!/bin/bash
#$ -cwd
#$ -N sqarr_full
#$ -j y
#$ -S /bin/bash
#$ -q highmem
#$ -l memoria_a_usar=460G
#$ -pe neworte 32
#$ -t 1-3
#$ -tc 1
#$ -V

# Is "only the two end loops couple" physical (twin-lead cancellation - the
# array runs down the centreline of a resonator hairpin fold) or a mesh
# artifact (all 46 junction ports sit in 3um gaps that a 4um fine mesh
# barely resolves)? Fine-mesh the WHOLE array uniformly and watch the
# middle loops. If they light up, the coupling is distributed and c is
# much larger; if they stay at ~1e-9, it's twin-lead cancellation.
#
# task: 1        2        3
#  array_1 min:  3.0      2.0      1.5   (um)   readout_1 fixed at 4um

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

AMIN_ARR=(3.0 2.0 1.5)
i=$((SGE_TASK_ID-1)); AMIN=${AMIN_ARR[$i]}
SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"; cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/fullarray/a${AMIN}"; mkdir -p "$OUTDIR"

( while true; do printf '%s ' "$(date +%H:%M:%S)"; free -g|awk 'NR==2{printf "used=%sG avail=%sG\n",$3,$7}'; sleep 30; done ) > "$OUTDIR/ram.log" 2>&1 &
RAMPID=$!; trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== fullarray task $SGE_TASK_ID: array min=${AMIN}um, host=$(hostname), $(date) ==="
timeout 43200 python -u build_and_eigenmode.py \
    --L-j-nH 600 --outdir "$OUTDIR" --palace-bin "$PALACE_BIN" --num-cpus "$NSLOTS" \
    --number-of-freqs 3 --starting-freq 7.6e9 \
    --bg-min-size-um 8 --resonator-bg-min-size-um 4 \
    --array-bg-min-size-um "$AMIN" --array-taper-dist-min-um 12 \
    --fillet-resolution 10 \
    2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}; kill $RAMPID 2>/dev/null
echo "=== task $SGE_TASK_ID done exit=$STATUS peakRAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m}' $OUTDIR/ram.log)G $(date) ==="
exit $STATUS
