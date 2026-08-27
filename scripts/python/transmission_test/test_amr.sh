#!/bin/bash
#$ -cwd
#$ -N test_amr
#$ -j y
#$ -S /bin/bash
#$ -q copahue,highmem
#$ -l memoria_a_usar=64G
#$ -l h=!compute-4-6&!compute-4-7&!compute-4-8
#$ -pe neworte 16
#$ -V

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg

PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/transmission_test"
cd "$SCRIPT_DIR"

OUTDIR=$HOME/circuit_design/transmission_test/pruebas
mkdir -p $OUTDIR

echo "===Eigenvalue test ==="
echo "=== Running on host: $(hostname) ==="

python build_and_eigenmode.py \
    --ground-pos 3.9527 \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS \
    --starting-freq 14.5e9 \
    --target-freq-ghz 15.75 \
    --number-of-freqs 5 \
    --amr-max-its 1 \
    --strap-taper-dist-max-um 20