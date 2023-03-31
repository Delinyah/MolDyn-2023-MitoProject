#!/bin/bash

# Set the base directory
base_dir=~/Project/OMM/structures-inner-new

# Loop through all the subfolders
for folder in "$base_dir"/*; do
    if [ -d "$folder" ]; then
        echo "Processing folder: $folder"
               
        # Run the Python script with the specified files
        python3 fppanel.py "$folder"
        echo "Completed processing folder: $folder"
    fi
done

# this script loops through all folders and calls a python script that makes a panel of four .png's that are in each subfolder
