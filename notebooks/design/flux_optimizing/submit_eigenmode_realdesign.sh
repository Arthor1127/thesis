#!/bin/bash
#$ -cwd
#$ -N flux_eigen_real
#$ -j y
#$ -S /bin/bash
#$ -q highmem,copahue,be_copahue
#$ -l memoria_a_usar=180G
#$ -pe neworte 32
#$ -V

# Eigenmode corroboration on the REAL target design (squid_array_coupling's
# build_design.py, verbatim defaults: N=23, squid_gap=25um, meandered readout
# at RES_LENGTH_MM=3.74575mm -> their own converged f_r=7.9355 GHz), not one
# of flux_optimizing's own sweep corners. See run_eigenmode_realdesign.py.
#
# 180G: sized between the small flux_optimizing corner (68G measured, 16
# cores) and the biggest-loop corner (183G measured, 24 cores) -- this
# design's per-loop geometry is closer to the small end (squid_gap=25 vs
# 14-40 swept) but the meandered 3.75mm resonator (jogs+fillets) adds more
# mesh than flux_optimizing's plain hairpin. Under copahue's 188G node ceiling
# with a little margin; bump to 220G (forces highmem) if this undershoots.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg QISKIT_METAL_HEADLESS=1
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

cd "$HOME/circuit_design/flux_optimizing"
mkdir -p eigen/realdesign
echo "=== flux_eigen_real host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
timeout 72000 python -u run_eigenmode_realdesign.py \
    --outdir eigen/realdesign \
    --palace-bin "$PALACE_BIN" --num-cpus "$NSLOTS" \
    2>&1 | tee -a eigen/realdesign/run.log
STATUS=${PIPESTATUS[0]}
echo "=== flux_eigen_real exit=$STATUS $(date) ==="
exit $STATUS
