#!/bin/bash
#$ -cwd
#$ -N qw_boxdiag
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=32G
#$ -l h=!compute-4-6&!compute-4-7&!compute-4-8
#$ -pe neworte 16
#$ -V

# DIAGNOSTIC: identify the unexplained 12.331 GHz mode.
#
# The quarter-wave baseline returned [5.054, 12.331, 15.126] GHz.
# Modes 1 and 3 are f1 and 3*f1 to within 0.2% - textbook open-short
# ladder. Mode 2 at 2.44*f1 belongs to no odd-harmonic ladder and is
# unexplained.
#
# It is not academic: the sweep's shift-invert targets for task 7
# (12.624 GHz) and task 9 (12.319 GHz) sit within 0.3 GHz and 0.012 GHz
# of it respectively, so --target-freq-ghz would select THIS mode rather
# than the real segment mode at those positions, silently corrupting two
# of the ten sweep points.
#
# THE TEST: box/cavity and other domain-supported modes scale with the
# chip and airbox dimensions; resonator modes scale with the trace
# length. This run changes ONLY the chip footprint (3.0x1.5mm ->
# 3.4x1.7mm, a ~13% linear increase) and leaves the resonator length
# untouched at 5.92905mm.
#
# READING THE RESULT, comparing 'All modes found' against the reference
# [5.054, 12.331, 15.126]:
#   * 5.054 and 15.126 stay put, 12.331 moves DOWN by roughly the
#     inverse of the size change (~13% -> ~10.9 GHz)
#       -> it is a domain/box mode. Options: size the box so it sits
#          outside 4-16 GHz, or keep it and filter post-hoc using
#          resonance.json's all_modes_ghz, selecting the mode nearest
#          the two-segment prediction while ignoring known box modes.
#   * ALL THREE move together
#       -> the change perturbed the resonator too (ground-plane
#          proximity), so the test is inconclusive; retry with a bigger
#          box change and check f1 stayed put first.
#   * 12.331 does not move at all
#       -> not a box mode. Next suspects are a slotline mode on the CPW
#          or something localised at the ShortToGround/OpenToGround
#          terminations. Visualise it in ParaView (mode index 1 in the
#          eigenmode .pvd) to see where the field actually lives.
#
# --starting-freq is set to 12.0 GHz here, deliberately targeting the
# unexplained mode's neighbourhood rather than f1, so the solver
# definitely returns it and its neighbours.

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg

PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/batch_test"
cd "$SCRIPT_DIR"

OUTDIR=$HOME/sweep_qw/boxdiag
mkdir -p $OUTDIR

echo "=== BOX-MODE DIAGNOSTIC: chip 3.0x1.5mm -> 3.4x1.7mm, L unchanged ==="
echo "=== Reference (3.0x1.5mm): [5.054, 12.331, 15.126] GHz ==="
echo "=== Running on host: $(hostname) ==="

python build_and_eigenmode_qw.py \
    --ground-pos none \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --total-length-mm 5.92905 \
    --chip-size-x-mm 3.4 --chip-size-y-mm 1.7 \
    --starting-freq 12.0e9 \
    --number-of-freqs 6 \
    --bg-min-size-um 8 --bg-max-size-um 100 \
    --target-freq-ghz 12.331
