#!/bin/bash
#$ -cwd
#$ -N sqarr_final
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=180G
#$ -pe neworte 32
#$ -t 1-6
#$ -tc 4
#$ -V

# Final runs, now that 1938540 showed only the two END loops couple (the
# array body screens the interior):
#   tasks 1-3: n=23, converge the end-loop p by shrinking the targeted
#              coupling-zone box (2.0 -> 1.0 -> 0.7 um) - this also pins the
#              relative sign of c_0 vs c_22
#   tasks 4-6: n_squids trend at box=2.0 um (n = 16, 10, 5); the far-end
#              loop stays at x=-1290, the near-end loop recedes from the
#              resonator, so c should fall with n
# All L_J = 600 nH, non-AMR (memory-safe).
#
# task:  1     2     3     4     5     6
#  n  :  23    23    23    16    10    5
#  box:  2.0   1.0   0.7   2.0   2.0   2.0  (um)

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

N_ARR=(23 23 23 16 10 5)
BOX_ARR=(2.0 1.0 0.7 2.0 2.0 2.0)
i=$((SGE_TASK_ID-1))
NSQ=${N_ARR[$i]}; BOX=${BOX_ARR[$i]}

SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"; cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/final/n${NSQ}_box${BOX}"; mkdir -p "$OUTDIR"

( while true; do printf '%s ' "$(date +%H:%M:%S)"; free -g|awk 'NR==2{printf "used=%sG avail=%sG\n",$3,$7}'; sleep 30; done ) > "$OUTDIR/ram.log" 2>&1 &
RAMPID=$!; trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== final task $SGE_TASK_ID: n=${NSQ} box=${BOX}um L_J=600, host=$(hostname), $(date) ==="
timeout 21600 python -u build_and_eigenmode.py \
    --L-j-nH 600 --n-squids "$NSQ" --outdir "$OUTDIR" \
    --palace-bin "$PALACE_BIN" --num-cpus "$NSLOTS" \
    --number-of-freqs 3 --starting-freq 7.6e9 \
    --bg-min-size-um 8 --resonator-bg-min-size-um 4 \
    --array-bg-min-size-um 4 --array-taper-dist-min-um 12 \
    --coupling-rect-min-um "$BOX" --coupling-rect-margin-um 150 \
    --fillet-resolution 10 \
    2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}; kill $RAMPID 2>/dev/null
echo "=== task $SGE_TASK_ID done exit=$STATUS peakRAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m}' $OUTDIR/ram.log)G $(date) ==="
exit $STATUS
