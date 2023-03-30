#!/bin/bash

# Set the base directory
base_dir=~/Project/OMM/structures-outer-new

# Loop through all the subfolders
for folder in "$base_dir"/*; do
    if [ -d "$folder" ]; then
        echo "Processing folder: $folder"
        # Set the path to the md.gro and md.xtc files
        gro_file="$folder/md.gro"
        xtc_file="$folder/md.xtc"
        
        # Run the Python script with the specified files
        python3 lipid_analysis.py "$gro_file" "$xtc_file"
    fi
done
