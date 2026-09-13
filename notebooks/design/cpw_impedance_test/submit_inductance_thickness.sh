#!/bin/bash
#$ -N cpwind_t
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
# (was copahue -- saturated by the flux_sweep campaign; highmem has slack and
# already ran the capacitance jobs fine)
#$ -l memoria_a_usar=24G
#$ -pe neworte 8
#$ -tc 6
#$ -o cpwind_t_task_$TASK_ID.log
#$ -t 1-12

# -o is a FLAT path on purpose (see submit_capacitance.sh's comment): SGE
# opens it before the job body's `mkdir -p` runs.

# Film-thickness-correction check (inductance side): FULL3D at a few
# thicknesses per representative (width, gap) point -- see
# gen_thickness_params.py. -t must match that script's printed row count
# (default 12 = 2 points x 6 thicknesses).
#
# Resource sizing is a first guess, bumped above the main sweep's since
# FULL3D geometry is more expensive to mesh; pilot task 1 alone first.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/cpw_impedance_test"
cd "$SCRIPT_DIR"
mkdir -p out_ind_thickness

LINE=$(sed -n "${SGE_TASK_ID}p" params_thickness.txt)
read -r W G THICK TAG <<< "$LINE"

( while true; do printf '%s ' "$(date +%H:%M:%S)"; free -g | awk 'NR==2{printf "used=%sG avail=%sG\n",$3,$7}'; sleep 30; done ) > "out_ind_thickness/ram_${SGE_TASK_ID}.log" 2>&1 &
RAMPID=$!; trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== task $SGE_TASK_ID : $LINE -> TAG=$TAG  host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 7200 python -u run_inductance.py \
    --outdir out_ind_thickness \
    --palace-bin "$PALACE_BIN" \
    --cpw-width-um "$W" \
    --cpw-gap-um "$G" \
    --full3d --metal-thickness-nm "$THICK" \
    --tag "$TAG" \
    --num-cpus "$NSLOTS" \
    2>&1 | tee "out_ind_thickness/run_${SGE_TASK_ID}.log"
STATUS=${PIPESTATUS[0]}
kill $RAMPID 2>/dev/null
echo "=== task $SGE_TASK_ID done exit=$STATUS peakRAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m"G"}' out_ind_thickness/ram_${SGE_TASK_ID}.log) $(date) ==="
exit $STATUS
