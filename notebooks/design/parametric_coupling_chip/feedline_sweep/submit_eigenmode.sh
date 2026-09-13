#!/bin/bash
#$ -cwd
#$ -N qc_eigen
#$ -j y
#$ -S /bin/bash
#$ -q highmem,copahue
#$ -l memoria_a_usar=150G
#$ -pe neworte 24
#$ -o qc_eigen_task.log
#$ -V

# Eigenmode cross-check for the cs14_cl300_ref driven-S21 pilot point: that
# run showed TWO weak resonant features (7.72GHz and 7.86GHz -- see the
# phase-residual analysis in this session). NO ports are created here (see
# build_and_eigenmode.py's docstring) so LP1/LP2 are plain open-ended pads --
# this exposes both the hanger resonator's own mode AND any bare
# feedline-cavity mode (LP1<->LP2 with no port loading) that the driven sim's
# lumped ports might otherwise mask or distort.
#
# 150G/24-core is a re-guess: the geometry changed substantially since the last
# successful run (job 1938822, 39.8G total on the OLD 2mm-feedline design) --
# the feedline is now 12mm (CHIP_HALF_SPAN_UM=6000um, to detune the bare
# feedline-cavity mode from 8GHz), so TL1/TL2's 8um-fine mesh now runs 6x
# longer. EPR interfaces are gone and readout_1 uses gap-only meshing now
# (both cut memory vs 1938822), but the much bigger domain could still net
# increase it -- watch PALACE's own "Estimated peak" line in out.log.
#
# 24 cores (not 32): compute-4-18 (the only host in the highmem/neworte
# hostgroup) has been seen holding another user's 8-slot job, leaving only 24
# genuinely free -- see this session's job 1938820 troubleshooting.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

COUPLING_SPACE_UM="${COUPLING_SPACE_UM:-14}"
COUPLING_LENGTH_UM="${COUPLING_LENGTH_UM:-300}"
STARTING_FREQ="${STARTING_FREQ:-6e9}"
NUMBER_OF_FREQS="${NUMBER_OF_FREQS:-4}"
AMR_MAX_ITS="${AMR_MAX_ITS:-3}"
AMR_TOL="${AMR_TOL:-1e-2}"

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_chip"
cd "$SCRIPT_DIR"
mkdir -p out_eig

echo "=== qc_eigen cs=${COUPLING_SPACE_UM}um cl=${COUPLING_LENGTH_UM}um host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u feedline_sweep/build_and_eigenmode.py \
    --coupling-space-um "$COUPLING_SPACE_UM" --coupling-length-um "$COUPLING_LENGTH_UM" \
    --outdir out_eig --palace-bin "$PALACE_BIN" --num-cpus "$NSLOTS" \
    --starting-freq "$STARTING_FREQ" --number-of-freqs "$NUMBER_OF_FREQS" \
    --amr-max-its "$AMR_MAX_ITS" --amr-tol "$AMR_TOL" \
    2>&1
STATUS=$?
echo "=== qc_eigen exit=$STATUS $(date) ==="
exit $STATUS
