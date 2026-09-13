#!/bin/bash
#$ -N cplrerun
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
#$ -l memoria_a_usar=90G
#$ -l h_rss=70G
#$ -l h_rt=9000
#$ -pe neworte 32
#$ -tc 1
#$ -o coupler_sweep/rerun_$TASK_ID.log
set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
export PALACE_BIN="$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace"
export OUTDIR="coupler_sweep"
cd ~/circuit_design/transmon_coupler
mkdir -p coupler_sweep/results

# coupler-geometry -> coupling-capacitance sweep. Each point: 2 CircTransmons +
# 1 CapNInterdigitalRound, pads held 60 um of bare neck apart (direct pad-pad
# ~0.5 fF, per pilot). Palace AMR (4 cycles, 4M-DoF cap) -- faster and better
# converged than a uniform fine mesh, and safely below the ~9M-DoF size where
# Palace's error estimator stalls. Aitken-extrapolate C_c per point afterwards.
LINE=$(sed -n "${SGE_TASK_ID}p" params_rerun.txt)
read -r FC FL FG FW FTG WC MODE <<< "$LINE"
export FINGER_COUNT="$FC" FINGER_LENGTH_UM="$FL" FINGER_GAP_UM="$FG"
export FINGER_WIDTH_UM="$FW" FINGER_TIP_GAP_UM="$FTG" WITH_COUPLER="$WC"
export SOLVER_ORDER=2 BG_MAX_UM=80 TRANSMON_RES_UM=14
if [ "$MODE" = "amr" ]; then
    # tamed AMR: gentle Dorfler + low DoF cap so even a crossing-cycle overshoot
    # stays ~3-4M (< the ~9M size where Palace's error estimator stalls, and
    # << RAM). iters high so the DoF cap, not the iter count, is the limiter.
    export AMR=1 AMR_ITERS=9 AMR_MAX_DOFS=3000000 AMR_TOL=1e-3 AMR_DORFLER=0.35 COUPLER_RES_UM=4.0
    SUF="_amr"
else
    export AMR=0 COUPLER_RES_UM=1.2
    SUF="_u1.2"
fi
export TAG="N${FC}_L${FL}_g${FG}_w${FW}_t${FTG}$([ "$WC" = 0 ] && echo _nocpl)${SUF}"

( while true; do printf "%s %s\n" "$(date +%s)" \
   "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
   sleep 20; done ) > "coupler_sweep/ram_rerun_${SGE_TASK_ID}.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

echo "=== rerun $SGE_TASK_ID  $TAG  $(date) ==="
python -u run_coupler_pair_cap.py 2>&1
echo "=== rerun $SGE_TASK_ID exit=$?  peak RSS: $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' coupler_sweep/ram_rerun_${SGE_TASK_ID}.log)  $(date) ==="
