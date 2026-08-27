#!/bin/bash

# Configuration
REMOTE_USER="arturo.ruiz"
REMOTE_HOST="10.73.25.223"
REMOTE_DIR="sweep_driven"
LOCAL_DIR="/home/ruiz/Documents/thesis/scripts/python/transmission_test"

# Ensure the local target directory exists
mkdir -p "$LOCAL_DIR"

echo "Locating s_sweep.npz files on $REMOTE_HOST..."

# Locate NPZ files on the cluster and process each match
ssh "${REMOTE_USER}@${REMOTE_HOST}" "find ~/${REMOTE_DIR} -type f -name 's_sweep.npz'" | while read -r filepath; do
    
    # Extract only the first 'pos_X.XXXXX' match from the path
    pos_folder=$(echo "$filepath" | grep -oE 'pos_[0-9]+\.[0-9]+' | head -n 1)
    
    if [[ -n "$pos_folder" ]]; then
        local_filename="${LOCAL_DIR}/${pos_folder}_s_sweep.npz"
        echo "Downloading ${pos_folder}_s_sweep.npz..."
        scp -q "${REMOTE_USER}@${REMOTE_HOST}:${filepath}" "$local_filename"
    fi
    
done

echo "=== All NPZ files successfully downloaded to $LOCAL_DIR ==="