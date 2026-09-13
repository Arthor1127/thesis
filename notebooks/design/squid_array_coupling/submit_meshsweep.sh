#!/bin/bash
#$ -cwd
#$ -N sqarr_msweep
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=160G
#$ -pe neworte 32
#$ -t 1-6
#$ -tc 4
#$ -V

# Non-AMR mesh-convergence sweep. AMR on this design blew past the 503G
# highmem node (46 lumped ports + AMR's retained meshes / estimation /
# rebalancing overhead ~ 3x a plain solve). A single non-AMR eigen solve
# is far lighter (a 3um-array solve was 37G / 0.75M elems), so convergence
# is checked the classic way instead: solve at 3 fixed fine-mesh densities
# and watch the resonator frequency AND the coupling c (c ~ sqrt(p), so p
# convergence is the real check - EPR_C_EXTRACTION.md Sec 6.2).
#
# task:  1        2        3        4        5        6
#  L_J:  600      600      600      400      400      400
#  readout/array fine-mesh min size (um):
#         5/4      3.5/3    2.5/2    5/4      3.5/3    2.5/2

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

LJ_ARR=(600 600 600 400 400 400)
RMIN_ARR=(5.0 3.5 2.5 5.0 3.5 2.5)
AMIN_ARR=(4.0 3.0 2.0 4.0 3.0 2.0)
i=$((SGE_TASK_ID-1))
LJ=${LJ_ARR[$i]}; RMIN=${RMIN_ARR[$i]}; AMIN=${AMIN_ARR[$i]}

SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"
cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/meshsweep/Lj${LJ}_r${RMIN}_a${AMIN}"
mkdir -p "$OUTDIR"

( while true; do printf '%s ' "$(date +%H:%M:%S)"; free -g|awk 'NR==2{printf "used=%sG avail=%sG\n",$3,$7}'; sleep 30; done ) > "$OUTDIR/ram.log" 2>&1 &
RAMPID=$!; trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== msweep task $SGE_TASK_ID: L_J=${LJ} readout=${RMIN}um array=${AMIN}um, host=$(hostname), $(date) ==="
timeout 21600 python -u build_and_eigenmode.py \
    --L-j-nH "$LJ" \
    --outdir "$OUTDIR" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    --number-of-freqs 2 \
    --starting-freq 6.5e9 \
    --bg-min-size-um 8 \
    --resonator-bg-min-size-um "$RMIN" \
    --array-bg-min-size-um "$AMIN" \
    --resonator-taper-dist-min-um 50 \
    --array-taper-dist-min-um 25 \
    --fillet-resolution 10 \
    2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}
kill $RAMPID 2>/dev/null
echo "=== task $SGE_TASK_ID done, exit=$STATUS, peak RAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m}' $OUTDIR/ram.log)G, $(date) ==="
exit $STATUS
