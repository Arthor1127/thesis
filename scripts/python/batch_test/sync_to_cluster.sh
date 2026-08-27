#!/bin/bash
# Sync .py and .sh files from the current local directory to the cluster.
#
# Run from your LOCAL machine, in the directory holding the scripts.
#
#   ./sync_to_cluster.sh <target_dir> [--go] [--delete]
#
# <target_dir> is relative to /home/arturo.ruiz on the cluster.
#
#   ./sync_to_cluster.sh circuit_design/batch_test          # DRY RUN
#   ./sync_to_cluster.sh circuit_design/batch_test --go     # actually transfer
#   ./sync_to_cluster.sh circuit_design/batch_test --go --delete
#                                                 # also remove cluster-side .py/.sh
#                                                 # that no longer exist locally
#
# Defaults to a DRY RUN deliberately. This project has had three separate
# incidents of a script silently not being updated on the cluster (a stale
# stageB_driven.sh ran the wrong solver order and mesh settings for a full
# job; build_and_eigenmode.py once lacked a flag the submit script passed,
# killing the job instantly). Seeing the file list before it moves is
# cheap insurance.

set -euo pipefail

CLUSTER_USER="${CLUSTER_USER:-arturo.ruiz}"
CLUSTER_HOST="${CLUSTER_HOST:-10.73.25.223}"
DRY="--dry-run"
DELETE=""
TARGET=""
for arg in "$@"; do
    case "$arg" in
        --go)     DRY="" ;;
        --delete) DELETE="--delete" ;;
        -h|--help)
            sed -n '2,24p' "$0"; exit 0 ;;
        -*)
            echo "unknown option: $arg" >&2; exit 1 ;;
        *)
            if [ -n "$TARGET" ]; then
                echo "target directory given twice: '$TARGET' and '$arg'" >&2
                exit 1
            fi
            TARGET="$arg"
            ;;
    esac
done

if [ -z "$TARGET" ]; then
    echo "ERROR: no target directory given." >&2
    echo >&2
    echo "  usage: $0 <target_dir> [--go] [--delete]" >&2
    echo "  e.g.:  $0 circuit_design/batch_test" >&2
    echo >&2
    echo "  <target_dir> is relative to /home/${CLUSTER_USER} on the cluster." >&2
    echo "  Deliberately has no default: this script overwrites files, and a" >&2
    echo "  silent default is how you sync into the wrong directory." >&2
    exit 1
fi

# Strip any leading/trailing slashes so both 'foo/bar' and '/foo/bar/'
# resolve to the same remote path rather than /foo/bar or a doubled slash.
TARGET="${TARGET#/}"
TARGET="${TARGET%/}"

DEST="${CLUSTER_USER}@${CLUSTER_HOST}:/home/${CLUSTER_USER}/${TARGET}/"

echo "=== local  : $(pwd)"
echo "=== remote : ${DEST}"
if [ -n "$DRY" ]; then
    echo "=== MODE   : DRY RUN (add --go to transfer)"
else
    echo "=== MODE   : LIVE TRANSFER"
fi
if [ -n "$DELETE" ]; then
    echo "=== DELETE : cluster-side .py/.sh with no local counterpart WILL be removed"
fi
echo

# FLAT sync - top-level .py and .sh only, no recursion. Omitting
# --include='*/' means --exclude='*' also blocks directories, so rsync
# never descends. This deliberately matches the verification block below,
# which checks only the top level; a recursive sync with a flat check
# would transfer files it never verified.
#
# -c compares by CHECKSUM rather than size+mtime. Slower, but size+mtime
# can miss an edit that leaves the byte count unchanged - exactly the
# failure mode that has caused stale-file incidents here. At these file
# sizes the cost is irrelevant.
rsync -avz -c $DRY $DELETE \
    --include='*.py' \
    --include='*.sh' \
    --exclude='*' \
    ./ "$DEST"

if [ -n "$DRY" ]; then
    echo
    echo "Dry run only - nothing was transferred. Re-run with --go."
    exit 0
fi

echo
echo "=== Verifying: comparing local vs remote checksums ==="
# Trust nothing. Compare md5 of every synced file on both ends and report
# mismatches explicitly, rather than assuming rsync's exit code is enough.
LOCAL_SUMS=$(find . -maxdepth 1 \( -name '*.py' -o -name '*.sh' \) -printf '%f\n' \
             | sort | while read -r f; do echo "$(md5sum < "$f" | cut -d' ' -f1)  $f"; done)

REMOTE_SUMS=$(ssh "${CLUSTER_USER}@${CLUSTER_HOST}" \
    "cd /home/${CLUSTER_USER}/${TARGET} 2>/dev/null && \
     for f in *.py *.sh; do [ -e \"\$f\" ] && echo \"\$(md5sum < \"\$f\" | cut -d' ' -f1)  \$f\"; done | sort")

MISMATCH=0
while read -r sum name; do
    remote_line=$(echo "$REMOTE_SUMS" | grep " ${name}\$" || true)
    if [ -z "$remote_line" ]; then
        echo "  MISSING on cluster : $name"
        MISMATCH=1
    else
        remote_sum=$(echo "$remote_line" | cut -d' ' -f1)
        if [ "$sum" != "$remote_sum" ]; then
            echo "  DIFFERS            : $name"
            MISMATCH=1
        fi
    fi
done <<< "$LOCAL_SUMS"

if [ "$MISMATCH" -eq 0 ]; then
    echo "  All local .py/.sh files match the cluster copies."
else
    echo
    echo "  ^^ Fix the above before submitting any job. If rsync keeps failing"
    echo "     on a specific file, fall back to writing it directly on the"
    echo "     cluster with a heredoc (cat > file << 'EOF' ... EOF), which has"
    echo "     been the reliable method in this project."
    exit 1
fi

echo
echo "NOTE: meandered_grounded.py is NOT synced by this script even if it"
echo "sits in this directory - it belongs to the qiskit-metal INSTALL path"
echo "(~/.conda/envs/qcg-quantum-design/lib/python3.11/site-packages/"
echo "qiskit_metal/qlibrary/tlines/), not the working directory. Copying it"
echo "here would have no effect; that exact mix-up already cost one round of"
echo "confusing overlap-sweep results. Update it on the cluster with:"
echo "  python -c \"import qiskit_metal.qlibrary.tlines.meandered_grounded as m; print(m.__file__)\""
echo "and copy over the path it prints."
