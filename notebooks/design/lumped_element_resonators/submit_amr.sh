#!/bin/bash
#$ -cwd
#$ -N lumped_amr
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=96G
#$ -pe neworte 16
#$ -V

# 2-iteration AMR test for the lumped-element LC resonator (interdigital cap
# + meander inductor). The base (uniform, no-AMR) mesh already ran locally:
# 456439 elements / 577660 global unknowns, 11.9G peak. The first AMR
# refinement pass pushed global unknowns to ~1.98M (H1+ND+RT), and that's
# where the local machine OOM'd mid-solve - moving to the cluster to finish
# that pass and one more, so the two iterations' indicator fields can be
# compared (same approach used for the squid_resonator_epr AMR study).
#
# 96G requested on copahue (not highmem): local base-mesh peak was 11.9G for
# 577k unknowns: iteration 1 more than triples that to ~1.98M unknowns, so
# if the setup/preconditioner memory scales roughly with unknown count,
# iteration 1 could plausibly need ~40G and iteration 2 (a further Dorfler
# refinement) somewhat more - 96G gives headroom without requesting anywhere
# near a whole node. copahue has multiple nodes in the 100-240G range
# (checked via `qstat -F memoria_a_usar` before submitting), so this doesn't
# require monopolizing the single highmem node.

set -uo pipefail

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
PALACE_BIN="$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace"

SCRIPT_DIR="$HOME/circuit_design/lumped_element_resonators"
cd "$SCRIPT_DIR"

echo "=== Running on host: $(hostname) ==="
echo "=== NSLOTS: ${NSLOTS:-unset} ==="
echo "=== Start time: $(date) ==="

OUTDIR="$SCRIPT_DIR/amr_run1"
mkdir -p "$OUTDIR"

timeout 7200 python -u build_and_eigenmode.py \
    --outdir "$OUTDIR" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus 16 \
    --number-of-freqs 3 \
    --amr-max-its 2 \
    2>&1 | tee "$OUTDIR/amr.log"
STATUS=${PIPESTATUS[0]}

echo "=== End time: $(date) ==="
if [ "$STATUS" -eq 124 ]; then
    echo "!!! TIMED OUT after 2 hours"
elif [ "$STATUS" -ne 0 ]; then
    echo "!!! FAILED with exit code $STATUS"
else
    echo "=== SUCCESS ==="
fi
