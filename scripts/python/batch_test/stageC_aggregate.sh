#!/bin/bash
#$ -cwd
#$ -N stageC
#$ -j y
#$ -S /bin/bash
#$ -q copahue
#$ -l memoria_a_usar=4G
#$ -pe neworte 1
#$ -hold_jid stageA,stageB
#$ -V

# Waits for ALL tasks of both stageA and stageB array jobs to finish
# (successfully or not) before running, via -hold_jid above.

conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH

export MPLBACKEND=Agg

cd "$HOME/circuit_design/batch_test"

echo "=== Stage C: aggregating results and generating plots ==="
python aggregate_and_plot.py --sweep-dir $HOME/sweep_v2

echo ""
echo "=== Full sweep complete ==="
echo "Plots: $HOME/sweep_v2/plots/"
echo "  baseline_S21.png"
echo "  S21_evolution_colormap.png"
echo "  S21_evolution_data.npz  (raw grid for further analysis)"
