#!/bin/bash
#$ -cwd
#$ -N test_driven8
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=64G
#$ -pe neworte 16
#$ -V

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg

python build_and_driven.py --ground-pos 5.906 --solver-order 2 \
    --outdir ~/sweep_v2_test/driven_5.906_v8 --palace-bin $PALACE_BIN --num-cpus 16
