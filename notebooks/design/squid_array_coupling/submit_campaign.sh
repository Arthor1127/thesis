#!/bin/bash
#$ -cwd
#$ -N sqarr_camp
#$ -j y
#$ -S /bin/bash
#$ -q highmem
#$ -l memoria_a_usar=450G
#$ -pe neworte 32
#$ -t 1-7
#$ -tc 1
#$ -V

# c-extraction campaign. Same AMR settings as the (L_J=600, n=23)
# convergence run in submit_converge.sh. Serial (-tc 1) on the single
# highmem node so peak RAM never overlaps between tasks.
#
#   L_J convergence : (400,23) vs the separate (600,23) conv run
#   n_squids trend  : n = 16, 10, 5 at each L_J, array anchored at
#                     pos_x=-1290um (fewer loops -> the array's near end
#                     recedes from the resonator -> coupling should drop,
#                     i.e. grow with n)
#
# task:  1        2        3        4        5        6        7
#  L_J:  400      600      600      600      400      400      400
#  n  :  23       16       10       5        16       10       5

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
export PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

LJ_ARR=(400 600 600 600 400 400 400)
N_ARR=( 23  16  10   5  16  10   5)
i=$((SGE_TASK_ID-1))
export L_J_NH=${LJ_ARR[$i]}
export N_SQUIDS_OVR=${N_ARR[$i]}
export AMR_ITERS=7
export AMR_MAX_DOFS=10000000
export AMR_DORFLER=0.5
export STARTING_FREQ_HZ=6.5e9
export NUMBER_OF_FREQS=2
export RESONATOR_MIN_UM=4.0
export ARRAY_MIN_UM=3.5

SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"
cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/campaign/Lj${L_J_NH}_n${N_SQUIDS_OVR}"
mkdir -p "$OUTDIR"

( while true; do
    printf '%s ' "$(date +%H:%M:%S)"
    free -g | awk 'NR==2{printf "used=%sG free=%sG avail=%sG\n",$3,$4,$7}'
    sleep 30
  done ) > "$OUTDIR/ram.log" 2>&1 &
RAMPID=$!
trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== campaign task $SGE_TASK_ID: L_J=${L_J_NH}nH n=${N_SQUIDS_OVR}, host=$(hostname), $(date) ==="
timeout 172800 python -u run_amr.py "$OUTDIR" 2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}
kill $RAMPID 2>/dev/null
echo "=== task $SGE_TASK_ID done, exit=$STATUS, $(date) ==="
awk '{for(i=1;i<=NF;i++) if($i ~ /^used=/){gsub(/[used=G]/,"",$i); if($i+0>m) m=$i+0}} END{print "peak node RAM: "m" G"}' "$OUTDIR/ram.log"
exit $STATUS
