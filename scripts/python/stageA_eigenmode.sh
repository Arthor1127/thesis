#!/bin/bash
#$ -cwd
#$ -N stageA
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=64G
#$ -l h=!compute-4-7&!compute-4-8
#$ -pe neworte 16
#$ -t 1-10
#$ -V

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

# Positions generated from build_design.py's own sweep_positions_mm() so
# this never drifts out of sync with the design module. "none" (the
# ungrounded baseline) is always task 1.
POSITIONS=(none $(python3 -c "from build_design import sweep_positions_mm; print(' '.join(str(x) for x in sweep_positions_mm(n_points=9)))"))
# 1 baseline + 9 numeric positions = 10 tasks -> matches #$ -t 1-10 above.
# If you change n_points in sweep_positions_mm(), update the -t range to
# (n_points + 1).

POS=${POSITIONS[$((SGE_TASK_ID-1))]}

OUTDIR=$HOME/sweep_v2/geom_${POS}
mkdir -p $OUTDIR

echo "=== Stage A task $SGE_TASK_ID: ground_pos=$POS mm ==="
python build_and_eigenmode.py \
    --ground-pos $POS \
    --outdir $OUTDIR \
    --palace-bin "$PALACE_BIN" \
    --num-cpus $NSLOTS
