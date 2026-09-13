#!/bin/bash
#$ -cwd
#$ -N flux_eigen
#$ -j y
#$ -S /bin/bash
#$ -q highmem,copahue,be_copahue
#$ -l memoria_a_usar=160G
#$ -pe neworte 32
#$ -V

# One eigenmode corroboration point (see run_eigenmode.py). Params come in via
# qsub -v CLR=..,SEG=..,GAP=..,SEMI=..,AW=..,TAG=..  (TAG names the output dir
# under eigen/).
#
# 32 slots (the neworte PE's per-host ceiling on copahue's biggest nodes,
# compute-4-14/15/16/17, 188G) is this account's standing default for a
# single solve (cluster-resource-request-sizing) -- eigenmode runs are heavy
# enough that the user asked for it explicitly here too.
#
# 160G is sized for small/mid loop geometries on a 188G copahue node. For the
# BIGGEST loop corner (squid_gap~40, arm_semi_length~24, arm_width~8) Palace
# itself measured "Estimated peak per-node memory usage" = 183.2G at only
# 24 cores/5.3M unknowns (job 1938721, memoria_a_usar was only 96G there --
# ran fine but was genuinely under-requested, a real oversubscription risk on
# a shared node). For configs at or near that corner, override with
# `qsub -v ...,submit_eigenmode.sh -l memoria_a_usar=220G` -- this forces
# placement onto the single highmem node (503G; no copahue node has 220G),
# which is the correct choice here, not a violation of "prefer copahue" --
# that guidance is for requests copahue can actually satisfy.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg QISKIT_METAL_HEADLESS=1
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

CLR="${CLR:-40}"
SEG="${SEG:-6}"
GAP="${GAP:-14}"
SEMI="${SEMI:-8}"
AW="${AW:-3}"
TAG="${TAG:-small}"

cd "$HOME/circuit_design/flux_optimizing"
mkdir -p "eigen/$TAG"
echo "=== flux_eigen [$TAG] clr=$CLR seg=$SEG gap=$GAP semi=$SEMI aw=$AW  host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 72000 python -u run_eigenmode.py \
    --route-clearance-um "$CLR" --segment-length-um "$SEG" \
    --squid-gap-um "$GAP" --arm-semi-length-um "$SEMI" --arm-width-um "$AW" \
    --outdir "eigen/$TAG" \
    --palace-bin "$PALACE_BIN" --num-cpus "$NSLOTS" \
    2>&1 | tee -a "eigen/$TAG/run.log"
STATUS=${PIPESTATUS[0]}
echo "=== flux_eigen [$TAG] exit=$STATUS $(date) ==="
if [ "$STATUS" -eq 0 ]; then
    python -u extract_eigen_M.py "eigen/$TAG" 2>&1 | tee -a "eigen/$TAG/extract.log"
    STATUS=${PIPESTATUS[0]}
fi
exit $STATUS
