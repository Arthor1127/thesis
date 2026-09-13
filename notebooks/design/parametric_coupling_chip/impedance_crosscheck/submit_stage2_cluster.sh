#!/bin/bash
#$ -N impz_stage2
#$ -cwd -S /bin/bash -V -j y
#$ -q copahue,be_copahue
#$ -l memoria_a_usar=64G
#$ -pe neworte 32
#$ -o stage2_cluster.log

# Stage 2 of the CPW lumped-port impedance cross-check, cluster-only (see
# QC_STATUS.md "Update 2026-09-11" and local-machine-heavy-fem-jobs-unsafe
# memory -- the local run OOM-killed an unrelated process at Order 4, so this
# moved here entirely, mesh generation included).
#
# Same measurement as driven_sweep/check_line_impedance.py (which got ~40 ohm):
# the thesis' own Si CPW cross-section (w=20/s=11.5/h=525/eps_r=11.45), a
# 2-port through line, Order 2 and Order 4. stage2_analyze.py recovers Z0(f)
# exactly from (S11,S21) via the NRW-style uniform-line inversion, no fitting.

set -uo pipefail
conda deactivate 2>/dev/null || true
source $HOME/spack/share/spack/setup-env.sh
spack load openmpi@4.1.8 %gcc@12.5.0
conda activate qcg-quantum-design
unset LD_LIBRARY_PATH
export MPLBACKEND=Agg

HERE="$HOME/circuit_design/parametric_coupling_chip/impedance_crosscheck"
cd "$HERE"

JULIA_BIN="$HOME/julia-1.12.7/bin/julia"
JULIA_PROJECT="$HERE/julia_mesh_env"
PALACE_BIN="${PALACE_BIN:-$HOME/spack/opt/spack/linux-x86_64/palace-0.16.0-wrb255mwkg6qlt4vsd6w6ne6fx3wfpdz/bin/palace}"
NP="${NP:-32}"
LINE_LEN_UM="${LINE_LEN_UM:-4000}"
ORDERS="${ORDERS:-2 4}"
FREQ_STEP_GHZ="${FREQ_STEP_GHZ:-2.0}"
TRACE_UM="${TRACE_UM:-20.0}"
GAP_UM="${GAP_UM:-11.5}"
GROUND_UM="${GROUND_UM:-250.0}"
REFN="${REFN:-1}"
OUTDIR="${OUTDIR:-stage2_out}"

OUT="$HERE/$OUTDIR"
rm -rf "$OUT"; mkdir -p "$OUT"
MESH="$OUT/mesh.msh"

echo "=== impz_stage2 host=$(hostname) NSLOTS=$NSLOTS TRACE_UM=$TRACE_UM GAP_UM=$GAP_UM GROUND_UM=$GROUND_UM REFN=$REFN OUTDIR=$OUTDIR $(date) ==="
echo "### generating mesh ###"
"$JULIA_BIN" --project="$JULIA_PROJECT" "$HERE/gen_mesh.jl" "$MESH" "$REFN" "$LINE_LEN_UM" "$TRACE_UM" "$GAP_UM" "$GROUND_UM"
ls -la "$MESH"

for order in $ORDERS; do
  echo; echo "################## STAGE 2 : order $order ##################"
  CFG="$OUT/cfg_o${order}.json"
  python3 - "$HERE/si_cpw_line.json" "$CFG" "$order" "$MESH" "$FREQ_STEP_GHZ" <<'PY'
import json, sys
src, dst, order, mesh, step = sys.argv[1:6]
c = json.load(open(src))
c["Solver"]["Order"] = int(order)
c["Model"]["Mesh"] = mesh
c["Problem"]["Output"] = f"postpro_o{order}"
c["Solver"]["Driven"]["Samples"][0]["FreqStep"] = float(step)
json.dump(c, open(dst, "w"), indent=2)
print("wrote", dst)
PY
  ( cd "$OUT" && "$PALACE_BIN" -np "$NP" "$CFG" ) \
    || { echo "!! order $order failed -- see above"; continue; }
done

echo; echo "results:"; find "$OUT" -name 'port-S.csv' -printf '  %p\n'
echo "=== impz_stage2 exit=$? $(date) ==="
