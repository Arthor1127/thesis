#!/bin/bash
#$ -N chaingeo
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
#$ -l memoria_a_usar=8G
#$ -l h_rss=6G
#$ -l h_rt=900
#$ -pe neworte 2
#$ -o chain_cap/geomcheck.log
set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
cd ~/circuit_design/transmon_coupler
python -u -c "
import matplotlib; matplotlib.use('Agg', force=True)
import build_chain
from build_design import _design_bounds_mm
d = build_chain.build_chain(n_transmons=6); d.rebuild()
print('components:', list(d.components.keys()))
print('junctions :', len(d.qgeometry.tables['junction']))
b = _design_bounds_mm(d)
print(f'bbox mm: X=[{b[0]:.3f},{b[2]:.3f}] Y=[{b[1]:.3f},{b[3]:.3f}]')
from SQDMetal.PALACE.Capacitance_Simulation import PALACE_Capacitance_Simulation
print('SQDMetal import OK')
print('GEOMCHECK OK')
"
