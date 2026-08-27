#!/bin/bash
#$ -cwd
#$ -N stageB
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=64G
#$ -l h=!compute-4-7&!compute-4-8
#$ -pe neworte 16
#$ -t 1-10
#$ -V

# Each task runs ONE continuous 4.5-20.5 GHz sweep (covering all four
# harmonics at once) for one grounding position - no per-harmonic
# splitting anymore, so 10 positions = 10 tasks.
#
# Does not depend on Stage A (fixed frequency grid, not per-point
# resonance), so this can run concurrently with Stage A.

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

POSITIONS=(none $(python3 -c "from build_design import sweep_positions_mm; print(' '.join(str(x) for x in sweep_positions_mm(n_points=9)))"))
# 1 baseline + 9 numeric positions = 10 tasks -> matches #$ -t 1-10 above.

POS=${POSITIONS[$((SGE_TASK_ID-1))]}

OUTDIR=$HOME/sweep_v2/driven_${POS}
mkdir -p $OUTDIR

echo "=== Stage B task $SGE_TASK_ID: ground_pos=$POS mm, continuous 4.5-20.5 GHz sweep ==="
python build_and_driven.py \
    --ground-pos $POS \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS
