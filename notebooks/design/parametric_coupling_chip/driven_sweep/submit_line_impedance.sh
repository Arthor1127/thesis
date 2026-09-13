#!/bin/bash
#$ -N qc_lineZ
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue
#$ -l memoria_a_usar=32G
#$ -pe neworte 32
#$ -o lineZ_$JOB_ID.log

# Bare-CPW-line impedance check, isolating where the cs14_cl300 design's
# Zin=20.9+5.5j Ohm (measured, away from resonance) actually comes from --
# the CoupledLineTee coupler, the LaunchpadWirebond tapers, or the plain CPW
# cross-section itself. See check_line_impedance.py's docstring.
#
# No resonator in either geometry -> no resonance to find -> a coarse 9-point
# frequency list (6-10GHz) is enough; much cheaper than the resonator sweep
# (32G/16-core is generous headroom, not a tight estimate).
#
# MODE is set by the caller (qsub -v MODE=with_launchpads|line_only), NOT
# read from a params file -- this is a single point per mode, not an array.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

: "${MODE:?MODE must be set via qsub -v MODE=with_launchpads (or line_only)}"

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_chip"
cd "$SCRIPT_DIR"
export OUTDIR="${OUTDIR:-out_line}"
export PORT_IMPEDANCE_R="${PORT_IMPEDANCE_R:-50}"
export TAG="${TAG:-$MODE}"
mkdir -p "$OUTDIR"

echo "=== qc_lineZ MODE=$MODE TAG=$TAG PORT_IMPEDANCE_R=$PORT_IMPEDANCE_R host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u driven_sweep/check_line_impedance.py 2>&1
ST=$?
echo "=== qc_lineZ MODE=$MODE exit=$ST $(date) ==="
exit $ST
