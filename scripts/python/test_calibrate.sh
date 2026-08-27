#!/bin/bash
#$ -cwd
#$ -N test_cal
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=16G
#$ -pe neworte 16
#$ -l h=!compute-4-7&!compute-4-8&!compute-4-9
#$ -V

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg

python build_and_eigenmode.py --ground-pos none --total-length-mm 14.9669 \
    --outdir ~/sweep_v2/geom_none_cal2 --palace-bin $PALACE_BIN --num-cpus 16
