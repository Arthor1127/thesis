#!/bin/bash
#$ -cwd
#$ -N sqarr_mswp2
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=170G
#$ -pe neworte 32
#$ -t 1-6
#$ -tc 4
#$ -V

# Mesh sweep 2: the c-relevant convergence variable (user 2026-08-29) is
# how well the 3um junction gap AND the vacuum around it are resolved -
# so vary the fine-mesh TAPER DISTANCE (how far the fine zone reaches into
# the surrounding air) at a fixed, moderate min size, rather than only
# shrinking min size. Also: nf 2 -> 3 and starting_freq 6.5 -> 7.6 GHz
# (mesh sweep 1's nf=2 / 6.5GHz gave "zero converged modes" on the finer
# meshes - the shift was too far below the ~8 GHz resonator).
#
# task:         1     2     3     4     5     6
#  L_J:         600   600   600   400   400   400
#  taper_min/max (um) for BOTH readout_1 and array_1:
#               30/200 80/400 160/700  30/200 80/400 160/700
# fixed: readout min 4um, array min 3um

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

LJ_ARR=(600 600 600 400 400 400)
TMIN_ARR=(30 80 160 30 80 160)
TMAX_ARR=(200 400 700 200 400 700)
i=$((SGE_TASK_ID-1))
LJ=${LJ_ARR[$i]}; TMIN=${TMIN_ARR[$i]}; TMAX=${TMAX_ARR[$i]}

SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"
cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/meshsweep2/Lj${LJ}_taper${TMIN}-${TMAX}"
mkdir -p "$OUTDIR"

( while true; do printf '%s ' "$(date +%H:%M:%S)"; free -g|awk 'NR==2{printf "used=%sG avail=%sG\n",$3,$7}'; sleep 30; done ) > "$OUTDIR/ram.log" 2>&1 &
RAMPID=$!; trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== mswp2 task $SGE_TASK_ID: L_J=${LJ} taper=${TMIN}/${TMAX}um, host=$(hostname), $(date) ==="
timeout 21600 python -u build_and_eigenmode.py \
    --L-j-nH "$LJ" \
    --outdir "$OUTDIR" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    --number-of-freqs 3 \
    --starting-freq 7.6e9 \
    --bg-min-size-um 8 \
    --resonator-bg-min-size-um 4 \
    --array-bg-min-size-um 3 \
    --resonator-taper-dist-min-um "$TMIN" \
    --resonator-taper-dist-max-um "$TMAX" \
    --array-taper-dist-min-um "$TMIN" \
    --array-taper-dist-max-um "$TMAX" \
    --fillet-resolution 10 \
    2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}
kill $RAMPID 2>/dev/null
echo "=== task $SGE_TASK_ID done, exit=$STATUS, peak RAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m}' $OUTDIR/ram.log)G, $(date) ==="
exit $STATUS
