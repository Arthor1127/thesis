#!/bin/bash
#$ -N qc_zoom
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue,highmem,be_copahue
#$ -l memoria_a_usar=40G
#$ -pe neworte 16
#$ -tc 4
#$ -o zoom_task_$TASK_ID.log
#$ -t 1-1

# Phase D: narrow-window driven S21 + complex hanger-notch fit, for the design
# points the eigenmode method can't resolve (see ../QC_STATUS.md). One point
# per params row: CS  CL  TAG  FR_SEED_GHZ  (FR_SEED from the eigenmode result's
# fr_ghz -- build the file with gen_zoom_params.py). PROM adaptive sampling over
# a ~20 MHz window centred on FR_SEED, complex notch fit (fit_qc_notch.py).
#
# Override the window with -v HALFWIN_MHZ=... ,NPTS=... ; the resonance
# linewidth is fr/Ql, so at Qc~1e5 and a dielectric-limited Qi~9e4 the FWHM is
# ~90 kHz -- keep the step (2*HALFWIN_MHZ*1e3/(NPTS-1) kHz) well below that.
# Default 10 MHz / 2001 pts = 10 kHz step.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export QT_QPA_PLATFORM=offscreen MPLBACKEND=Agg
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"

SCRIPT_DIR="$HOME/circuit_design/parametric_coupling_chip"
cd "$SCRIPT_DIR"

sleep $(( (SGE_TASK_ID - 1) * 8 ))

PARAMS_FILE="${PARAMS_FILE:-driven_sweep/params_qc_zoom.txt}"
LINE=$(sed -n "${SGE_TASK_ID}p" "$PARAMS_FILE")
read -r CS CL TAG FRSEED <<< "$LINE"
export COUPLING_SPACE_UM="$CS" COUPLING_LENGTH_UM="$CL" TAG PALACE_BIN
export FR_SEED_GHZ="$FRSEED"
export HALFWIN_MHZ="${HALFWIN_MHZ:-10}" NPTS="${NPTS:-2001}"
export ADAPTIVE_MAX_SAMPLES="${ADAPTIVE_MAX_SAMPLES:-400}"
export SOLVER_ORDER="${SOLVER_ORDER:-1}" SOLVER_TOL="${SOLVER_TOL:-1e-8}" SOLVER_MAXITS="${SOLVER_MAXITS:-300}"
export MESH_FACTOR="${MESH_FACTOR:-0.75}"
export OUTDIR="${OUTDIR:-out_qc_zoom}"
mkdir -p "$OUTDIR"

( while true; do printf "%s %s\n" "$(date +%s)" \
    "$(ps -o rss= --ppid $$ -p $$ 2>/dev/null | awk '{s+=$1}END{print s/1024}')"; \
    sleep 5; done ) > "$OUTDIR/ram_${SGE_TASK_ID}.log" 2>&1 &
MON=$!; trap "kill $MON 2>/dev/null" EXIT

RUN_LOG="$OUTDIR/run_${SGE_TASK_ID}.log"
echo "=== qc_zoom task $SGE_TASK_ID : cs=$CS cl=$CL fr=$FRSEED -> TAG=$TAG host=$(hostname) NSLOTS=$NSLOTS $(date) ==="
python -u driven_sweep/run_qc_driven_zoom.py 2>&1 | tee "$RUN_LOG"
ST=${PIPESTATUS[0]}
kill $MON 2>/dev/null
echo "=== qc_zoom task $SGE_TASK_ID ($TAG) exit=$ST  palace: $(grep -iE 'estimated peak.*node' "$RUN_LOG" | tail -1)  $(date) ==="
exit $ST
