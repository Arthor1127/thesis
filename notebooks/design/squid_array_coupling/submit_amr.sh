#!/bin/bash
#$ -cwd
#$ -N sqarr_amr
#$ -j y
#$ -S /bin/bash
#$ -q highmem
#$ -l memoria_a_usar=200G
#$ -pe neworte 32
#$ -V

# AMR run (2 iterations) on the SQUID-array + readout-resonator design.
#
# L_J is raised to 600 nH so the whole SQUID-array phonon band drops well
# below the readout resonator (band top ~ 1/sqrt(L_J): the converged L_J=15nH
# run put it near 30 GHz, so 600 nH -> ~4.7 GHz, vs the ~8 GHz resonator).
# The eigen window therefore STARTS ABOVE the phonon band (6.5 GHz) so the
# solver targets the resonator directly and AMR refines for it, not for the
# array. identify_resonator_mode then picks it by sum|p_j| << 1.
#
# highmem (compute-4-18, ~503G, idle), 200G, 32 ranks. maxvmem on the
# earlier n23 AMR was 112G, so 200G is comfortable headroom.

set -uo pipefail

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
export PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"
export AMR_ITERS=2
export L_J_NH=600
export STARTING_FREQ_HZ=6.5e9
export NUMBER_OF_FREQS=4

SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"
cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/amr_Lj600"
mkdir -p "$OUTDIR"

echo "=== SQUID-array AMR: host=$(hostname), NSLOTS=$NSLOTS, L_J=${L_J_NH}nH, AMR_ITERS=$AMR_ITERS, $(date) ==="

timeout 21600 python -u run_amr.py "$OUTDIR" 2>&1 | tee "$OUTDIR/run.log"
STATUS=${PIPESTATUS[0]}
echo "=== AMR done, exit=$STATUS, $(date) ==="
if [ "$STATUS" -eq 124 ]; then
    echo "!!! TIMED OUT (6h) - partial per-iteration data is under $OUTDIR/*/iteration*/"
fi
exit $STATUS
