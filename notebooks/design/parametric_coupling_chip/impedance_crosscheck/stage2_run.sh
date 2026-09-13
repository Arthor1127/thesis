#!/usr/bin/env bash
# STAGE 2 -- the thesis' own CPW cross-section (Si, w=20 um, s=11.5 um), a
# 2-port THROUGH line (one 2-fin lumped port per end). Generates the mesh with
# Gmsh.jl, then solves the driven sweep at Order 2 and Order 4. Recovers Z0 and
# gamma EXACTLY per frequency from (S11, S21) -- stage2_analyze.py. Same
# measurement as ../driven_sweep/check_line_impedance.py (which got ~40 ohm).
#
#   PALACE_BIN=/home/ruiz/repo/spack/.../bin/palace \
#   JULIA_PROJECT=/home/ruiz/PalaceForCQED \
#   IMPEDANCE_LOCAL_OK=1 ./stage2_run.sh
#
# Knobs: REFINEMENT (default 1), TRACE_LEN_UM (4000), ORDERS ("2 4"),
#        FREQ_STEP_GHZ (2.0 -> 13 solves; drop to 0.5 once it looks sane).

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/_guard.sh"

JULIA_BIN="${JULIA_BIN:-julia}"
JULIA_PROJECT="${JULIA_PROJECT:-/home/ruiz/PalaceForCQED}"   # any env with Gmsh.jl
REFINEMENT="${REFINEMENT:-1}"
LINE_LEN_UM="${LINE_LEN_UM:-4000}"
ORDERS="${ORDERS:-2 4}"
FREQ_STEP_GHZ="${FREQ_STEP_GHZ:-2.0}"

OUT="$HERE/stage2_out"
rm -rf "$OUT"; mkdir -p "$OUT"
MESH="$OUT/mesh.msh"

echo "### generating mesh (Gmsh.jl via project $JULIA_PROJECT) ###"
"$JULIA_BIN" --project="$JULIA_PROJECT" "$HERE/gen_mesh.jl" "$MESH" "$REFINEMENT" "$LINE_LEN_UM"
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
  ( cd "$OUT" && MEMLOG="$OUT/mem_o${order}.log" palace_solve "$CFG" ) \
    || { echo "!! order $order failed/killed -- see $OUT/mem_o${order}.log"; continue; }
done

echo; echo "results:"; find "$OUT" -name 'port-S.csv' -printf '  %p\n'
echo "next: python $HERE/stage2_analyze.py"
