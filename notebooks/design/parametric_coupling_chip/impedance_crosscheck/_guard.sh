# Shared safety preamble for the LOCAL impedance cross-check runs.
# Sourced by stage1_run_palace_cpw.sh and stage2_run.sh.
#
# The thesis standing rule is: never run FEM/PALACE solves on this desktop
# (g15, 14 GB RAM + ~4 GB swap) -- a PALACE solve once drove it into swap
# exhaustion, OOM-killed the encrypted-vault daemon, and forced a hard reboot.
# These runs are a deliberate, user-authorised exception for a DELIBERATELY
# TINY problem (single short CPW line, coarse mesh). The guards below make a
# blow-up fail cleanly (std::bad_alloc) instead of taking the machine down.

set -euo pipefail

if [ "${IMPEDANCE_LOCAL_OK:-}" != "1" ]; then
  echo "REFUSING TO RUN: set IMPEDANCE_LOCAL_OK=1 to confirm you want a local PALACE solve." >&2
  echo "  (standing rule: FEM solves go to the CNEA cluster; this is a one-off exception)" >&2
  exit 3
fi

# Hard virtual-memory cap PER RANK -> malloc failure, not OOM-kill.
# User authorised all 8 cores; memory is the only cap. This tiny problem
# partitions to well under 1 GiB/rank, so 5 GiB/rank * 8 = generous ceiling
# that still can't collectively exhaust a 14 GB machine before the mem-logger
# + manual watch catch it (sustained swap growth = kill now).
VMEM_CAP_KB="${VMEM_CAP_KB:-5242880}"      # 5 * 1024 * 1024, per rank
ulimit -v "$VMEM_CAP_KB"
ulimit -c 0                                # no core dumps

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
NP="${NP:-8}"
RUN_TIMEOUT="${RUN_TIMEOUT:-2400}"          # 40 min wall cap per solve

PALACE_BIN="${PALACE_BIN:?set PALACE_BIN to the local palace wrapper}"

echo "=== guard: vmem cap $((VMEM_CAP_KB/1024)) MiB | np=$NP | omp=$OMP_NUM_THREADS | timeout ${RUN_TIMEOUT}s ==="
echo "--- memory before ---"; free -m
# background memory logger; the caller traps EXIT to kill it
_MEMLOG="${MEMLOG:-mem_$(date +%s).log}"
( while true; do date +%s | tr -d '\n'; echo -n ' '; free -m | awk '/Mem:/{print "used="$3" avail="$7} /Swap:/{print "swap_used="$3}' | paste -sd' '; sleep 3; done ) > "$_MEMLOG" 2>&1 &
_MEMPID=$!
trap 'kill $_MEMPID 2>/dev/null || true' EXIT
echo "memory log -> $_MEMLOG (watch it: sustained swap growth = Ctrl-C NOW)"

palace_solve() {  # palace_solve <config.json> [extra palace args...]
  local cfg="$1"; shift
  echo ">>> timeout ${RUN_TIMEOUT} $PALACE_BIN -np $NP $cfg $*"
  timeout "$RUN_TIMEOUT" "$PALACE_BIN" -np "$NP" "$cfg" "$@"
}
