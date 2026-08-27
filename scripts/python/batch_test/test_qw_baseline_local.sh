#!/bin/bash
# Quarter-wave BASELINE test - local machine.
#
# Same role as the half-wave design's own baseline check: confirms
# whether TOTAL_LENGTH_MM=5.92905mm in build_design_quarterwave.py
# actually gives f1=5GHz for this real geometry (fillets, lead-in
# straights, and terminations all add second-order effects the ideal
# formula misses - the half-wave design needed two correction passes
# before converging).
#
# Uses build_and_eigenmode_qw.py, which is build_and_eigenmode.py with
# its only change being the import - it points at
# build_design_quarterwave.build_design instead of the half-wave one.
# Same CLI, same validated fast settings (8/100um mesh - the measured
# sweet spot; NOT the 20um mesh that gave spurious modes because it
# could not resolve the CPW cross-section; no AMR, since SQDMetal's own
# 'mesh_refinement' key is UniformLevels not AMR and the real AMR flags
# --amr-max-its/--amr-tol are unvalidated).
#
# This design has no TL/launchpads, so it is far smaller than the
# half-wave sweep's eigenmode runs - expect well under the ~30min the
# half-wave baseline took at this same mesh setting.
#
# Usage: place next to build_design_quarterwave.py and
# build_and_eigenmode_qw.py, then just run it.
#   bash test_qw_baseline_local.sh

set -e

PALACE_BIN="${PALACE_BIN:-/home/ruiz/repo/spack/opt/spack/linux-zen3/palace-develop-blpbtxe43ne2or3yo7gou6ycmdwdxuzl/bin/palace}"
NUM_CPUS="${NUM_CPUS:-4}"

python build_and_eigenmode_qw.py \
    --ground-pos none \
    --outdir ~/sweep_qw_test_local \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NUM_CPUS" \
    --total-length-mm 5.92905 \
    --starting-freq 5.0e9 \
    --bg-min-size-um 8 --bg-max-size-um 100 \
    --target-freq-ghz 5.0

echo ""
echo "Check the BASELINE CHECK line above. If |%off| > 5%, the script"
echo "itself will suggest a corrected --total-length-mm - feed that back"
echo "into TOTAL_LENGTH_MM in build_design_quarterwave.py and rerun."
