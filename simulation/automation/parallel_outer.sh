#!/bin/bash

# Directory where PDB files are stored
pdb_dir="/martini/delinyah/Project/OMM/structures-outer-new/"

# Maximum number of parallel jobs
max_jobs=10

# Function to check the number of running jobs
function running_jobs {
    jobs | wc -l
}

# Iterate over all PDB files
for pdb_file in $pdb_dir/*.pdb
do
  # Remove the .pdb extension from the file name
  pdb_name=$(basename $pdb_file .pdb)

  # Run the simulation script in the background for each PDB file
  ./mainone.sh $pdb_name &

  # If the number of running jobs has reached the limit, wait for any job to finish
  while (( $(running_jobs) >= max_jobs ))
  do
    sleep 1
  done
done

# Wait for all background jobs to finish
wait
