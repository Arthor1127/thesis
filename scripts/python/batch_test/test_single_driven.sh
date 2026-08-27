#!/bin/bash
#$ -cwd
#$ -N test_driven9
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

# v9: first test of the split port/resonator background mesh (see
# build_and_driven.py). Ports (TL/LP1/LP2) keep the validated 15/150um;
# resonator1 goes independently coarser at 30/300um, since it no longer
# needs port-grade resolution along its whole ~15mm length now that the
# grounding strap gets its own separately-targeted fine mesh box. New
# outdir/job name (v9, not v8) so this doesn't overwrite the previous
# test_driven8 output - compare the two logs directly once this finishes.
python build_and_driven.py --ground-pos 5.906 --solver-order 2 \
    --outdir ~/sweep_v2_test/driven_5.906_v9 --palace-bin $PALACE_BIN --num-cpus 16 \
    --resonator-bg-min-size-um 30 --resonator-bg-max-size-um 300
