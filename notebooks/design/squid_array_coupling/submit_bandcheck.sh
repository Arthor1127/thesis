#!/bin/bash
#$ -cwd
#$ -N sqarr_band
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=64G
#$ -pe neworte 32
#$ -V

# Phonon-band check at L_J = 600 nH: solve the low modes on a mesh fine
# enough to actually resolve the 3 um junction gaps (array_bg_min = 3 um,
# not the 6-8 um of the L_J scan, which returned spurious hot-spot modes
# with numerically-zero port-EPR). No AMR. Confirms the whole array phonon
# band sits below the ~8 GHz readout resonator with a clear gap.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/squid_array_coupling"
cd "$SCRIPT_DIR"
OUTDIR="$SCRIPT_DIR/bandcheck_Lj600"
mkdir -p "$OUTDIR"

echo "=== band check L_J=600nH: host=$(hostname), NSLOTS=$NSLOTS, $(date) ==="

timeout 5400 python -u build_and_eigenmode.py \
    --L-j-nH 600 \
    --outdir "$OUTDIR" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NSLOTS" \
    --number-of-freqs 12 \
    --starting-freq 3.5e9 \
    --bg-min-size-um 6 \
    --resonator-bg-min-size-um 6 \
    --array-bg-min-size-um 3 \
    --fillet-resolution 10 \
    2>&1 | tee "$OUTDIR/run.log"
echo "=== band check done, exit=${PIPESTATUS[0]}, $(date) ==="
