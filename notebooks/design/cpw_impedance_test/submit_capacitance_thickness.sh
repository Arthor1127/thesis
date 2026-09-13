#!/bin/bash
#$ -N cpwcap_t
#$ -cwd -S /bin/bash -V -j y
#$ -q highmem
#$ -l memoria_a_usar=24G
#$ -l h_rss=18G
#$ -l h_rt=3600
#$ -pe neworte 8
#$ -tc 6
#$ -o cpwcap_t_task_$TASK_ID.log
#$ -t 1-12

# -o is a FLAT path on purpose (see submit_capacitance.sh's comment): SGE
# opens it before the job body's `mkdir -p` runs.

# Film-thickness-correction check (capacitance side): FULL3D (real metal
# thickness) at a few thicknesses per representative (width, gap) point -- see
# gen_thickness_params.py for the grid and DesignHelper.
# compute_cpw_geometric_parameters's Eqs. (B6)-(B7) for the analytic
# correction being tested. -t must match that script's printed row count
# (default 12 = 2 points x 6 thicknesses).
#
# Resource sizing is a first guess, bumped above the main (zero-thickness)
# sweep's since FULL3D geometry is more expensive to mesh; pilot task 1 alone
# first.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/cpw_impedance_test"
cd "$SCRIPT_DIR"
mkdir -p out_cap_thickness

ulimit -v $((32 * 1024 * 1024)) || true

LINE=$(sed -n "${SGE_TASK_ID}p" params_thickness.txt)
read -r W G THICK TAG <<< "$LINE"
export CPW_WIDTH_UM="$W" CPW_GAP_UM="$G" TAG PALACE_BIN
export FULL3D=1 METAL_THICKNESS_NM="$THICK"
export OUTDIR="out_cap_thickness"

( while true; do printf "%s %s\n" "$(date +%s)" \
    "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
    sleep 5; done ) > "out_cap_thickness/ram_${SGE_TASK_ID}.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

echo "=== task $SGE_TASK_ID : $LINE -> TAG=$TAG  host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u run_capacitance.py 2>&1
ST=$?
kill $MON 2>/dev/null
echo "=== task $SGE_TASK_ID exit=$ST  peak RSS: $(awk 'BEGIN{m=0}{if($2>m)m=$2}END{printf "%.0f MB",m}' out_cap_thickness/ram_${SGE_TASK_ID}.log 2>/dev/null) ($(date)) ==="
exit $ST
