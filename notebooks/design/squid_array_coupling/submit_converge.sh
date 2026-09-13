#!/bin/bash
#$ -cwd
#$ -N sqarr_conv
#$ -j y
#$ -S /bin/bash
#$ -q highmem
#$ -l memoria_a_usar=480G
#$ -pe neworte 32
#$ -V

# Convergence AMR for L_J=600 nH, n_squids=23.
#
# v2: the v1 run (Dorfler 0.6, nf=3, chip margin 4x) hit Palace's own
# memory estimate of 394G at 9.4M DoF and stopped after only 3 AMR
# iterations (freq step still 0.58%, not < 0.5%). Fixes:
#   - chip margin 4x -> 2.5x in build_design.py (~40% less mesh volume)
#   - number_of_freqs 3 -> 2
#   - Dorfler 0.6 -> 0.5, iterations 5 -> 7 (gentler, more steps)
#   - DoF cap 18M -> 10M  (Palace est ~42G/M-DoF -> ~420G, safe on 503G)
# Base mesh still seeds readout_1 @ 4um / array_1 @ 3.5um per the user's
# boundary-output inspection. Frequency tracked per iteration; a background
# logger watches node RAM.
#
# A node-RAM logger runs in the background -> ram.log.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
export PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

export L_J_NH=600
export N_SQUIDS_OVR=23
export AMR_ITERS=7
export AMR_MAX_DOFS=10000000
export AMR_DORFLER=0.5
export STARTING_FREQ_HZ=6.5e9
export NUMBER_OF_FREQS=2
export RESONATOR_MIN_UM=4.0
export ARRAY_MIN_UM=3.5

SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"
cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/conv_Lj600_n23_v2"
mkdir -p "$OUTDIR"

# --- background node-RAM logger ---
( while true; do
    printf '%s ' "$(date +%H:%M:%S)"
    free -g | awk 'NR==2{printf "used=%sG free=%sG avail=%sG\n",$3,$4,$7}'
    top -bn1 | grep -E "palace-x86_64|python" | awk '{print "   "$12" RSS="$6/1048576"G"}' | head -4
    sleep 20
  done ) > "$OUTDIR/ram.log" 2>&1 &
RAMPID=$!
trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== conv AMR L_J=600 n=23: host=$(hostname), NSLOTS=$NSLOTS, $(date) ==="
free -g | head -2

timeout 172800 python -u run_amr.py "$OUTDIR" 2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}
kill $RAMPID 2>/dev/null

echo "=== conv AMR done, exit=$STATUS, $(date) ==="
echo "--- peak node RAM seen ---"
awk '{for(i=1;i<=NF;i++) if($i ~ /^used=/){gsub(/[used=G]/,"",$i); if($i+0>m) m=$i+0}} END{print m" G"}' "$OUTDIR/ram.log"
exit $STATUS
