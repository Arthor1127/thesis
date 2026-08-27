#!/bin/bash
# Local (laptop) runner for the inductance smoke test.
#
# Nothing like the cluster's setup_env.sh - no spack, no proxy, no
# offscreen Qt workarounds beyond the two exports below. Just conda plus
# a palace binary.
#
#   ./run_local_test.sh                 # mesh only, no solve
#   ./run_local_test.sh --solve         # mesh + solve
#   ./run_local_test.sh --solve --gui   # also open the GMSH GUI

set -euo pipefail

CONDA_ENV="${CONDA_ENV:-qcg-quantum-design}"
OUTDIR="${OUTDIR:-$HOME/palace_ind_test}"
NUM_CPUS="${NUM_CPUS:-4}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

EXTRA_ARGS=("--mesh-only")
for arg in "$@"; do
    case "$arg" in
        --solve) EXTRA_ARGS=() ;;
        --gui)   EXTRA_ARGS+=("--view-gmsh") ;;
        *)       EXTRA_ARGS+=("$arg") ;;
    esac
done

# Headless rendering. qiskit_metal force-switches the matplotlib backend at
# import time, which is why the Python script calls matplotlib.use("Agg")
# again AFTER importing SQDMetal rather than relying on this alone.
export QT_QPA_PLATFORM=offscreen
export MPLBACKEND=Agg

echo "=================================================="
echo " Local inductance smoke test"
echo "=================================================="

# --- conda ---
if [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV" ]]; then
    CONDA_BASE="$(conda info --base 2>/dev/null)" || {
        echo "ERROR: conda not on PATH."; exit 1; }
    # shellcheck disable=SC1091
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV"
fi
echo "conda env: ${CONDA_DEFAULT_ENV:-<none>}"
echo "python:    $(which python)"

# --- palace binary ---
if [[ -z "${PALACE_BIN:-}" ]]; then
    PALACE_BIN="$(command -v palace || true)"
fi
if [[ -z "$PALACE_BIN" || ! -x "$PALACE_BIN" ]]; then
    echo ""
    echo "!!! No palace binary found."
    echo "    Set PALACE_BIN to its path, e.g.:"
    echo "      export PALACE_BIN=\$HOME/spack/opt/.../bin/palace"
    echo ""
    echo "    You can still run the mesh-only check without it, which"
    echo "    verifies the geometry and the generated config JSON."
    if [[ " ${EXTRA_ARGS[*]:-} " != *"--mesh-only"* ]]; then
        echo "    Refusing to --solve without a solver. Exiting."
        exit 1
    fi
    PALACE_BIN="palace"
else
    echo "palace:    $PALACE_BIN"
    if command -v mpirun >/dev/null 2>&1; then
        echo "mpirun:    $(which mpirun)"
    else
        echo "mpirun:    NOT FOUND - Palace needs MPI to run in parallel."
    fi
fi

# --- SQDMetal ---
echo -n "SQDMetal:  "
python - <<'PY'
import importlib.util, os, sys
spec = importlib.util.find_spec("SQDMetal")
if spec is None:
    print("NOT IMPORTABLE"); sys.exit(1)
root = os.path.dirname(spec.origin)
print(root)
f = os.path.join(root, "PALACE", "Inductance_Simulation.py")
print("           Inductance_Simulation.py:",
      "present" if os.path.exists(f) else "MISSING - git pull your clone")
PY

mkdir -p "$OUTDIR"
echo "outdir:    $OUTDIR"
echo "--------------------------------------------------"

python "$SCRIPT_DIR/test_inductance_local.py" \
    --palace-bin "$PALACE_BIN" \
    --outdir "$OUTDIR" \
    --num-cpus "$NUM_CPUS" \
    "${EXTRA_ARGS[@]}"
