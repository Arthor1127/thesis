#!/bin/bash
#$ -cwd
#$ -N flux_sweep_hi
#$ -j y
#$ -S /bin/bash
#$ -q highmem,copahue,be_copahue
#$ -l memoria_a_usar=44G
#$ -pe neworte 12
#$ -t 1-150
#$ -tc 60
#$ -V

# highmem listed FIRST: compute-4-18 (503 G, 32 slots) is non-preemptible and
# idle at night, unlike be_copahue on the same node.
#
# memoria_a_usar = 32 G (NOT enforced on this cluster -- h_vmem is INFINITY -- so
# it is purely a placement hint; undersizing it lets a node oversubscribe and the
# kernel OOM-killer can then take out other users' jobs). Measured maxvmem is
# 8-12 G for small/mid geometries; the big-loop corner (gap 40, semi 24, aw 8)
# is ~5-6 M tets -> estimated ~16-20 G with MaxSize=200. 32 G covers the worst
# case with margin and still fits 4 jobs on a 188 G copahue node. cpuINgpu
# dropped -- its usable nodes are small/flaky and not worth the OOM risk.

# 5-D flux-threading sweep: 4*5*6*5*4 = 2400 grid points over
#   route_clearance_um, segment_length_um (inter-loop lead), squid_gap_um (loop
#   height), arm_semi_length_um (loop half-width), arm_width_um (lead width).
# Array metal IN, array near (elbow-side) end tied to ground. Each point is a
# magnetostatic Palace solve at the CONVERGED boundary-layer mesh (~3.4 M tets,
# <1% mesh error -- flat 3 um was -15%); no paraview; GMRES. shard_<K>.csv is
# flushed per point and resumes past finished points on resubmit.
#
# GAP-FILL pass: 150 shards x 16 pts, 12 cores + 48 G (headroom for any big-
# loop points that OOM-failed the main 8c/32G pass). Resumes past done points. Queues copahue + be_copahue (preemptible,
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

echo "=== flux_sweep shard $SGE_TASK_ID/150: host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 72000 python -u sweep_flux.py \
    --shard "$SGE_TASK_ID" --nshards 150 \
    --outdir "$OUT" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    2>&1 | tee -a "$OUT/run_$K.log"
STATUS=${PIPESTATUS[0]}
kill $RAMPID 2>/dev/null
echo "=== shard $SGE_TASK_ID exit=$STATUS peakRAM $(awk '{for(i=1;i<=NF;i++)if($i~/used=/){gsub(/[used=G]/,"",$i);if($i+0>m)m=$i+0}}END{print m"G"}' $OUT/ram_$K.log) $(date) ==="
exit $STATUS
