#!/bin/bash
#$ -N cpwZ0chk
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem,copahue
#$ -l memoria_a_usar=16G
#$ -l h_rt=14400
#$ -pe neworte 8
#$ -tc 1
#$ -o z0chk_task_$TASK_ID.log
#$ -t 1-2

# First attempt (h_rt=3600, tc=2/concurrent) hit two problems: task 2 crashed
# on "./temp.sh: Text file busy" (both tasks run `cd` into the same
# $SCRIPT_DIR, and PALACE_Capacitance_Simulation.run() writes a wrapper
# script there by relative path -- concurrent tasks collide on it), and task 1
# got silently killed by the 1hr h_rt wall-clock limit mid-solve (no exit
# trailer -- this w20/g11.5 cross-section is bigger than this pipeline's
# previously-validated test points, needs more mesh/DOFs, and evidently more
# than an hour). Fixed here: tc=1 serializes the two tasks (no temp.sh
# collision) and h_rt raised to 4 hours.

# Independent Z0 cross-check for the parametric_coupling_chip design's
# cpw_width=20um/cpw_gap=11.5um cross-section, motivated by an unexplained
# driven-S-parameter result: Zin extracted via Z0_ref(1+S11)/(1-S11) on a bare
# CPW line (no launchpad, no coupler) converged (after ruling out mesh
# resolution, vacuum height, and lateral ground-plane width -- all tested and
# eliminated, see check_line_impedance.py's convergence/lateral sweeps) to
# ~40.5 Ohm, not the ~50.3 Ohm both DesignHelper and SQDMetal's own CpwParams
# analytic formula predict for the same cross-section (which also matches the
# 11.45 permittivity the actual PALACE config uses -- confirmed, not a
# material mismatch either).
#
# This runs a genuinely DIFFERENT extraction method: capacitance-per-length
# via PALACE_Capacitance_Simulation (the same pipeline this directory already
# validated for other cross-sections -- w8_g5 gave +9.4% vs analytic, a normal
# finite-length end-fringing residual, not a 20%+ gap). Two runs, same
# geometry, only the dielectric differs (real silicon vs vacuum,
# DIELECTRIC_MATERIAL=vacuum -> Material("") -> er=1 in SQDMetal). For a
# quasi-TEM line L' is purely geometric, so this sidesteps the SEPARATE,
# still-unresolved ~40,000x bug in this repo's own inductance sim entirely:
#
#   L' = 1/(c0^2 * C'_vacuum)          (geometric, dielectric-independent)
#   Z0 = sqrt(L'/C'_real) = 1/(c0*sqrt(C'_vacuum*C'_real))
#
# If this Z0 lands near 50 Ohm, the FEM/mesh physics is fine and the ~40.5 Ohm
# driven-sim number is an S-parameter EXTRACTION artifact (candidate: SQDMetal
# builds CPW LumpedPorts as two separate "fin" elements per port, not a single
# lumped element -- worth checking against Palace's own S-parameter definition
# for that port type). If this ALSO lands near 40 Ohm, the discrepancy is
# something more fundamental shared by both methods.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/cpw_impedance_test"
cd "$SCRIPT_DIR"
mkdir -p out_z0check

LINE=$(sed -n "${SGE_TASK_ID}p" params_z0check.txt)
read -r W G DIEL TAG <<< "$LINE"
export CPW_WIDTH_UM="$W" CPW_GAP_UM="$G" DIELECTRIC_MATERIAL="$DIEL" TAG PALACE_BIN
export OUTDIR="out_z0check"

echo "=== task $SGE_TASK_ID : $LINE -> TAG=$TAG DIEL=$DIEL host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u run_capacitance.py 2>&1
ST=$?
echo "=== task $SGE_TASK_ID exit=$ST ($(date)) ==="
exit $ST
