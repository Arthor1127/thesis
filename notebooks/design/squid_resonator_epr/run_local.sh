#!/bin/bash
# Run build_and_eigenmode.py locally (no cluster/SGE submission).
#
# Usage:
#   ./run_local.sh --outdir ~/sweep_v2/run1 [--L-transmon-nH 10] [--L-probe-nH 15] \
#       [--res-length-mm 3.74575] [--num-cpus 16] [any other build_and_eigenmode.py flag...]
#
# Any flag not recognized here is passed straight through to
# build_and_eigenmode.py, so e.g. --junction-mesh-min-um, --number-of-freqs,
# --amr-max-its etc. all work unchanged.
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P)"
SETUP_ENV="$HOME/Documents/thesis/scripts/python/setup_env.sh"

if [[ ! -f "$SETUP_ENV" ]]; then
    echo "Error: could not find $SETUP_ENV" >&2
    exit 1
fi

# --- Environment (conda + spack MPI + Palace binary path) ---
source "$SETUP_ENV"

# --- Defaults ---
OUTDIR=""
NUM_CPUS=16
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --num-cpus)
            NUM_CPUS="$2"
            shift 2
            ;;
        *)
            EXTRA_ARGS+=("$1")
            shift
            ;;
    esac
done

if [[ -z "$OUTDIR" ]]; then
    echo "Error: --outdir is required" >&2
    exit 1
fi
mkdir -p "$OUTDIR"

echo "=== Running build_and_eigenmode.py locally ==="
echo "  outdir:      $OUTDIR"
echo "  num-cpus:    $NUM_CPUS"
echo "  palace-bin:  $PALACE_BIN"
echo "  extra args:  ${EXTRA_ARGS[*]}"
echo "==============================================="

python "$SCRIPT_DIR/build_and_eigenmode.py" \
    --outdir "$OUTDIR" \
    --palace-bin "$PALACE_BIN" \
    --num-cpus "$NUM_CPUS" \
    "${EXTRA_ARGS[@]}"
