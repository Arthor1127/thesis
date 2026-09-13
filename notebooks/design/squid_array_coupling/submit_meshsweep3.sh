#!/bin/bash
#$ -cwd
#$ -N sqarr_mswp3
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=180G
#$ -pe neworte 32
#$ -t 1-4
#$ -tc 4
#$ -V

# Mesh sweep 3: TARGETED fine mesh on the two loops that actually carry the
# coupling (array far end + near end), via fine_mesh_in_rectangle boxes with
# a vacuum margin - instead of a global fine taper over the whole 5 mm of
# metal perimeter, which blew the mesh up to 33 M unknowns in sweep 2.
# Sweep the box min size to converge p (hence c ~ sqrt(p)).
#
# task:              1     2     3     4
#  L_J:              600   600   400   400
#  coupling box min: 3.0   2.0   3.0   2.0  (um), boxes +/-150 um

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

LJ_ARR=(600 600 400 400)
BOX_ARR=(3.0 2.0 3.0 2.0)
i=$((SGE_TASK_ID-1))
LJ=${LJ_ARR[$i]}; BOX=${BOX_ARR[$i]}

SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"
cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/meshsweep3/Lj${LJ}_box${BOX}"
mkdir -p "$OUTDIR"

( while true; do printf '%s ' "$(date +%H:%M:%S)"; free -g|awk 'NR==2{printf "used=%sG avail=%sG\n",$3,$7}'; sleep 30; done ) > "$OUTDIR/ram.log" 2>&1 &
RAMPID=$!; trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== mswp3 task $SGE_TASK_ID: L_J=${LJ} box_min=${BOX}um, host=$(hostname), $(date) ==="
timeout 21600 python -u build_and_eigenmode.py \
    --L-j-nH "$LJ" \
    --outdir "$OUTDIR" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    --number-of-freqs 3 \
    --starting-freq 7.6e9 \
    --bg-min-size-um 8 \
    --resonator-bg-min-size-um 4 \
    --array-bg-min-size-um 4 \
    --resonator-taper-dist-min-um 40 \
    --array-taper-dist-min-um 15 \
    --coupling-rect-min-um "$BOX" \
    --coupling-rect-margin-um 150 \
    --fillet-resolution 10 \
    2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}
kill $RAMPID 2>/dev/null
echo "=== task $SGE_TASK_ID done, exit=$STATUS, peak RAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m}' $OUTDIR/ram.log)G, $(date) ==="
exit $STATUS
