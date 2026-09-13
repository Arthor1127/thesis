#!/bin/bash
#$ -cwd
#$ -N flux_sweep
#$ -j y
#$ -S /bin/bash
#$ -q highmem,copahue,be_copahue
#$ -l memoria_a_usar=32G
#$ -pe neworte 8
#$ -t 1-300
#$ -tc 60
#$ -V

# highmem listed FIRST: compute-4-18 (503 G, 32 slots) is non-preemptible and
# idle at night, unlike be_copahue on the same node.
#
# memoria_a_usar = 44 G, 12 cores. This resource is NOT enforced here
# (h_vmem = INFINITY) -- it is a placement hint, so undersizing it lets a node
# oversubscribe and the kernel OOM-killer can hit other users. Measured peak
# was ~31 G at 8c/MaxSize-200; ~40 G expected at 12c/MaxSize-300. 44 G fits on
# every 62 G+ node with margin. cpuINgpu dropped (small/flaky nodes).

# 5-D flux-threading sweep: 4*5*6*5*4 = 2400 grid points over
#   route_clearance_um, segment_length_um (inter-loop lead), squid_gap_um (loop
#   height), arm_semi_length_um (loop half-width), arm_width_um (lead width).
# Array metal IN, array near (elbow-side) end tied to ground. Each point is a
# magnetostatic Palace solve at the CONVERGED boundary-layer mesh (~3.4 M tets,
# <1% mesh error -- flat 3 um was -15%); no paraview; GMRES. shard_<K>.csv is
# flushed per point and resumes past finished points on resubmit.
#
# 300 shards x 8 points; 8 cores/job. Queues copahue + be_copahue (preemptible,
# fine -- shards resume) + cpuINgpu. ~5 min/point -> ~9-12 h wall depending on
# free slots.

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
OUT="sweep"
mkdir -p "$OUT"
K=$(printf '%03d' "$SGE_TASK_ID")

( while true; do printf '%s ' "$(date +%H:%M:%S)"; free -g | awk 'NR==2{printf "used=%sG avail=%sG\n",$3,$7}'; sleep 60; done ) > "$OUT/ram_$K.log" 2>&1 &
RAMPID=$!
# always tear down the scratch workdir, even on qdel/timeout (SIGTERM) -- an
# orphaned mid-solve mesh+field dir is 20-100 MB
trap 'kill $RAMPID 2>/dev/null; rm -rf "$SCRIPT_DIR/$OUT/work_$K"' EXIT

echo "=== flux_sweep shard $SGE_TASK_ID/300: host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 72000 python -u sweep_flux.py \
    --shard "$SGE_TASK_ID" --nshards 300 \
    --outdir "$OUT" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    2>&1 | tee -a "$OUT/run_$K.log"
STATUS=${PIPESTATUS[0]}
kill $RAMPID 2>/dev/null
echo "=== shard $SGE_TASK_ID exit=$STATUS peakRAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m"G"}' $OUT/ram_$K.log) $(date) ==="
exit $STATUS
