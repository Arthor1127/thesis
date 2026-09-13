#!/bin/bash
#$ -N qc_fdln
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem,copahue
#$ -l memoria_a_usar=90G
#$ -pe neworte 8
#$ -tc 4
#$ -o fdln_task_$TASK_ID.log
#$ -t 1-6

# Feedline-length sweep, eigenmode, at a FIXED coupling point (cs=14um,
# cl=300um). The question this answers:
#
#   Job 1938824 (12mm feedline) found 4 modes -- 6.883, 7.671, 8.291, 8.738 GHz
#   -- and the per-mode E-field plots showed only ONE of them (8.291GHz) is
#   resonator-localized; the other three are feedline/launchpad-dominated. So
#   lengthening the feedline from 2mm to 12mm did NOT push its modes away from
#   the resonator, it packed MORE of them around it. That is the expected
#   behaviour in hindsight: an open-open cavity's harmonic spacing is v/(2L),
#   i.e. INVERSELY proportional to length, so a longer feedline has a DENSER
#   mode ladder. The earlier "lengthen it to detune it" reasoning was backwards.
#
# This sweep measures the trend directly instead of arguing about it. Because
# feedline harmonics scale as 1/L while the resonator's own mode does not, a
# plot of every returned mode vs feedline length separates them with no field
# analysis at all: the ~flat line near 8.3GHz is the resonator, the fanning
# 1/L family is the feedline. Whichever length maximises the gap between the
# flat line and its nearest neighbour is the one to design to.
#
# CONFOUND, and the probe that resolves it: chip_size_x auto-grows with the
# feedline, so feedline length and simulation-box size move together and a mode
# that shifts across the sweep could be either a feedline harmonic OR a box /
# substrate mode. That distinction matters, because only the first is fixed by
# choosing a feedline length -- and the 12mm field plots showed strong field
# along the chip's top/bottom edges, which looks box-like. Worse, the simple
# v/(2L) harmonic model does NOT reproduce the measured 12mm modes (it predicts
# 5.03/10.06GHz; the feedline-dominated modes actually came back at 6.883,
# 7.671 and 8.738GHz), so the analytic picture is already known to be wrong.
# Task 2 (hs3000_bigbox) therefore repeats the 6mm-feedline point inside the
# SAME 14.8mm box as the 12mm point: compare it against task 4 (6mm feedline,
# 8.8mm box) and task 1 (12mm feedline, 14.8mm box). Modes that track the box
# move with task 1; modes that track the feedline move with task 4.
#
# Geometry note: for chip_half_span_um >= 3000 the resonator is byte-identical
# (same single fold, x in [150,1850]um, same depth), so those points differ
# ONLY in feedline length. hs=2000 is near-identical; hs=1500 is a single fold
# with slightly different leg split; hs=1000 folds twice (its snake x_max is
# capped to stay clear of the closer launchpads). Total resonator LENGTH is the
# same everywhere, so f_r should barely move -- but if the hs=1000 point looks
# like an outlier, that different fold pattern is the first thing to suspect.
#
# No AMR (it segfaulted a rank in job 1938821 and has not been re-tested since
# the EPR removal) and no EPR interfaces (not needed for mode frequencies).
#
# Resources: 8 cores x 3 concurrent so the whole array fits inside highmem's
# 24 free slots on compute-4-18 WITHOUT depending on copahue, which has been
# saturated/overloaded by other users all evening. 90G covers the largest
# point with margin (the 12mm/24-rank run peaked at 57.6G total; smaller chips
# scale down roughly with chip area). Longest task is listed FIRST so it starts
# in the first wave -- SGE dispatches array tasks in order.
#
# -o is a FLAT path on purpose (SGE opens it before the job body's mkdir runs).

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_chip"
cd "$SCRIPT_DIR"

LINE=$(sed -n "${SGE_TASK_ID}p" feedline_sweep/params_feedline.txt)
read -r HALF_SPAN CHIP_X TAG <<< "$LINE"
OUTDIR="out_fdln/$TAG"
mkdir -p "$OUTDIR"

# CHIP_X of "auto" lets build_design size the chip to the launchpads; an
# explicit value pins the box size instead (the decoupling probe, see header).
CHIP_X_ARG=""
if [ "$CHIP_X" != "auto" ]; then
    CHIP_X_ARG="--chip-size-x-mm $CHIP_X"
fi

echo "=== qc_fdln task $SGE_TASK_ID ($TAG): half_span=${HALF_SPAN}um (feedline $(echo "$HALF_SPAN" | awk '{print $1*2/1000}')mm) chip_x=${CHIP_X} host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u feedline_sweep/build_and_eigenmode.py \
    --coupling-space-um 14 --coupling-length-um 300 \
    --chip-half-span-um "$HALF_SPAN" $CHIP_X_ARG \
    --outdir "$OUTDIR" --palace-bin "$PALACE_BIN" --num-cpus "$NSLOTS" \
    --starting-freq 5e9 --number-of-freqs 6 \
    --amr-max-its 0 \
    2>&1
ST=$?
echo "=== qc_fdln task $SGE_TASK_ID ($TAG) exit=$ST $(date) ==="
exit $ST
