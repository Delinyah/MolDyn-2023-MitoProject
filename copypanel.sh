#!/bin/bash

# Set the base directory
base_dir=~/Project/OMM/structures-inner-new

# Set the destination directory for panel images
dest_dir=~/Project/OMM/panels-gf

# Create the destination directory if it doesn't exist
mkdir -p "$dest_dir"

# Loop through all the subfolders
for folder in "$base_dir"/*; do
    if [ -d "$folder" ]; then
        echo "Processing folder: $folder"
               
        # Run the Python script with the specified files
        python3 fppanel.py "$folder"
        echo "Completed processing folder: $folder"

        # Copy the panel image to the destination directory
        basename=$(basename "$folder")
        cp "$folder/$basename-panel.png" "$dest_dir/"
    fi
done
