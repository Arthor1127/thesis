#!/bin/bash
#$ -cwd
#$ -N flux_arr
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=32G
#$ -pe neworte 32
#$ -t 1-2
#$ -tc 2
#$ -V

# Magnetostatic flux-per-amp of the readout hairpin through the 23 SQUID-loop
# apertures. Two tasks, identical except for the array metal:
#   task 1 : array metal IN   (reproduces inductance_sim.ipynb: ~11.8 Phi0/mA,
#                               GMRES, ~570 iters -- the array's open loops are
#                               topological handles that stall AMS)
#   task 2 : array metal OUT  (loop apertures are vacuum, ground opened by a
#                               RectCutout)
# If task 2's flux_per_amp is within a few % of task 1, the sim can drop the
# array metal and run ~5x faster on plain CG. Otherwise the array stays in.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/flux_optimizing"
cd "$SCRIPT_DIR"

if [ "$SGE_TASK_ID" -eq 1 ]; then
    OUT="cluster_flux_array";   FLAG=""
else
    OUT="cluster_flux_noarray"; FLAG="--no-array"
fi
mkdir -p "$OUT"

( while true; do printf '%s ' "$(date +%H:%M:%S)"; free -g | awk 'NR==2{printf "used=%sG avail=%sG\n",$3,$7}'; sleep 30; done ) > "$OUT/ram.log" 2>&1 &
RAMPID=$!; trap 'kill $RAMPID 2>/dev/null' EXIT

echo "=== flux task $SGE_TASK_ID ($OUT): host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 3600 python -u run_flux.py \
    --outdir "$OUT" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    $FLAG \
    2>&1 | tee "$OUT/run.log"
STATUS=${PIPESTATUS[0]}
kill $RAMPID 2>/dev/null
echo "=== task $SGE_TASK_ID done exit=$STATUS peakRAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m"G"}' $OUT/ram.log) $(date) ==="
exit $STATUS
