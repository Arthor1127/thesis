#!/bin/bash
#$ -N cplpilot
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
#$ -l memoria_a_usar=80G
#$ -l h_rss=64G
#$ -l h_rt=14400
#$ -pe neworte 32
#$ -tc 1
#$ -o coupler_sweep/pilot_$TASK_ID.log
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

LINE=$(sed -n "${SGE_TASK_ID}p" params_pilot.txt)
read -r FC FL FG FW FTG WC MODE <<< "$LINE"
export FINGER_COUNT="$FC" FINGER_LENGTH_UM="$FL" FINGER_GAP_UM="$FG"
export FINGER_WIDTH_UM="$FW" FINGER_TIP_GAP_UM="$FTG" WITH_COUPLER="$WC"
export SOLVER_ORDER=2 BG_MAX_UM=80 TRANSMON_RES_UM=14
if [ "$MODE" = "amr" ]; then
    export AMR=1 AMR_ITERS=4 AMR_MAX_DOFS=4000000 AMR_TOL=5e-3 COUPLER_RES_UM=2.5
    export TAG="N${FC}_L${FL}_g${FG}_w${FW}_t${FTG}$([ "$WC" = 0 ] && echo _nocpl)_amr"
else
    export AMR=0 COUPLER_RES_UM=1.0
    export TAG="N${FC}_L${FL}_g${FG}_w${FW}_t${FTG}$([ "$WC" = 0 ] && echo _nocpl)_u1.0"
fi

( while true; do printf "%s %s\n" "$(date +%s)" \
   "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
   sleep 15; done ) > "coupler_sweep/ram_pilot_${SGE_TASK_ID}.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

echo "=== pilot $SGE_TASK_ID  $TAG  $(date) ==="
python -u run_coupler_pair_cap.py 2>&1
echo "=== pilot $SGE_TASK_ID exit=$?  peak RSS: $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' coupler_sweep/ram_pilot_${SGE_TASK_ID}.log) $(date) ==="
