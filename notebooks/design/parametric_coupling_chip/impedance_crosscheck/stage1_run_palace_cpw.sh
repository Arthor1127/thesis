#!/usr/bin/env bash
# STAGE 1 -- sanity: does THIS palace binary reproduce its own documented CPW
# lumped-port benchmark (Z0 = 56.02 ohm, examples/cpw)? Uses the pre-built mesh
# shipped with the Palace source -- NO new geometry, NO gmsh. Runs the stock
# Order-2 config and an Order-4 copy (the docs say higher order is needed to
# get the lumped-port CPW reflection down).
#
#   PALACE_BIN=/home/ruiz/repo/spack/.../bin/palace \
#   IMPEDANCE_LOCAL_OK=1 ./stage1_run_palace_cpw.sh
#
# Output: ./stage1_out/{order2,order4}/postpro/lumped_uniform/port-S.csv
# Then:   python stage1_analyze.py

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/_guard.sh"

PALACE_SRC="${PALACE_SRC:-/home/ruiz/repo/palace-src}"
CPW_DIR="$PALACE_SRC/examples/cpw"
[ -f "$CPW_DIR/cpw_lumped_uniform.json" ] || { echo "no $CPW_DIR/cpw_lumped_uniform.json" >&2; exit 2; }
[ -f "$CPW_DIR/mesh/cpw_lumped_0.msh" ]   || { echo "pre-built mesh missing" >&2; exit 2; }

OUT="$HERE/stage1_out"
rm -rf "$OUT"; mkdir -p "$OUT"

# work from a copy so the palace-src tree stays clean
cp -r "$CPW_DIR" "$OUT/cpw_src"
cd "$OUT/cpw_src"

python3 - <<'PY'
import json, re, pathlib
base = json.loads(re.sub(r'//.*', '', pathlib.Path('cpw_lumped_uniform.json').read_text()))
# trim the sweep to keep it cheap: 3 GHz steps, no extra 17 GHz point-save,
# no field saves. We only need |S11|/arg(S11) on a coarse grid to back out Zin.
base['Solver']['Driven'] = {
    "Samples": [{"Type": "Linear", "MinFreq": 3.0, "MaxFreq": 30.0, "FreqStep": 3.0}],
}
base['Domains']['Postprocessing'] = {}
for order, sub in ((2, 'order2'), (4, 'order4')):
    c = json.loads(json.dumps(base))
    c['Solver']['Order'] = order
    c['Problem']['Output'] = f"postpro_{sub}"
    pathlib.Path(f"cfg_{sub}.json").write_text(json.dumps(c, indent=2))
    print("wrote", f"cfg_{sub}.json")
PY

for sub in order2 order4; do
  echo; echo "################## STAGE 1 : $sub ##################"
  MEMLOG="$OUT/mem_$sub.log" palace_solve "cfg_$sub.json" || { echo "!! $sub solve failed/killed"; continue; }
done

echo; echo "done. results:"
find "$OUT" -name 'port-S.csv' -printf '  %p\n'
echo "next: python $HERE/stage1_analyze.py"
