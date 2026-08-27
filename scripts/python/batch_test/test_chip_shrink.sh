#!/bin/bash
#$ -cwd
#$ -N eig_shrink
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=64G
#$ -l h=!compute-4-6&!compute-4-7&!compute-4-8
#$ -pe neworte 16
#$ -V

# Validates the CHIP SHRINK in isolation. Everything except the chip
# footprint matches the known-good 8/100um reference run, so any change
# in f1 is attributable to the shrink alone.
#
# REFERENCE (8/100um, chip 4.8 x 2.4mm, L=11.9907mm):
#     f1 = 4.9169 GHz,  2.8M unknowns,  1952s
#
# THIS RUN: same mesh, same length, chip Y shrunk 2.4 -> 1.8mm with the
# centre moved -1.0 -> -0.7mm so the content stays centred. That keeps
# 250um of ground below otg2 (~40 gap-widths, far beyond where coplanar
# fields are screened) while cutting chip volume to ~75%.
#
# WHAT TO COMPARE:
#   f1 shift   - if <0.1% the shrink is free and should be adopted in
#                build_design.py's CHIP_SIZE_Y_MM / CHIP_CENTER_Y_MM
#                defaults. A larger shift means the ground plane or
#                airbox was doing real work and 1.8mm is too tight; try
#                2.0 / -0.8 instead.
#   unknowns   - expect roughly 75% of 2.8M
#   Preconditioner line in the timing report - the number that matters,
#                since it is ~75% of runtime and scales worse than
#                linearly with unknowns.
#
# Note the length here is 11.9907mm, matching the REFERENCE run, NOT the
# 11.7914mm from the 6/80um run. Comparing against the right baseline
# matters more than using the latest length.

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

OUTDIR=$HOME/sweep_v2/eig_shrink_test
mkdir -p $OUTDIR

echo "=== CHIP SHRINK validation: chip Y 2.4 -> 1.8mm, centre -1.0 -> -0.7 ==="
echo "=== Reference to beat: f1=4.9169 GHz, 2.8M unknowns, 1952s ==="
echo "=== Running on host: $(hostname) ==="

python build_and_eigenmode.py \
    --ground-pos none \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --total-length-mm 11.9907 \
    --starting-freq 4.9e9 \
    --bg-min-size-um 8 --bg-max-size-um 100 \
    --chip-size-y-mm 1.8 --chip-center-y-mm -0.7 \
    --target-freq-ghz 5.0
