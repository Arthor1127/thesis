#!/bin/bash

# Configuration
REMOTE_USER="arturo.ruiz"
REMOTE_HOST="10.73.25.223"
REMOTE_DIR="sweep_driven"
LOCAL_DIR="$HOME/Downloads/transmission_test"

# Ensure the local directory exists
mkdir -p "$LOCAL_DIR"

echo "Locating s_params.png files on $REMOTE_HOST..."

# ssh into the cluster, find all matching images, and read them line by line
ssh "${REMOTE_USER}@${REMOTE_HOST}" "find ~/${REMOTE_DIR} -type f -name 's_params.png'" | while read -r filepath; do
    
    # Extract the position folder name (e.g., 'pos_3.77077') from the remote path
    pos_folder=$(echo "$filepath" | grep -o 'pos_[0-9.]*')
    
    if [[ -n "$pos_folder" ]]; then
        # Construct the new local file name
        local_filename="${LOCAL_DIR}/${pos_folder}_s_params.png"
        
        echo "Downloading data for ${pos_folder}..."
        # Download the file silently (-q)
        scp -q "${REMOTE_USER}@${REMOTE_HOST}:${filepath}" "$local_filename"
    fi
    
done

echo "=== All sweeps successfully downloaded to $LOCAL_DIR ==="