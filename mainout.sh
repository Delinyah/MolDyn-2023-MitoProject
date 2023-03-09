#!/bin/bash

# Help function 
function help {
  echo "Automated Workflow for Outer Mitochondrial Membrane Simulation Modeling"
  echo "Title: mainone (for MArtinize-INsanify Outer-membraNE)
  echo "Author: Delinyah C. Koning (as of March 9, 2023)"
  echo ""
  echo "Usage: omain.sh <pdb_code> [-nt <number_of_threads>]"
  echo "Don't forget to give permission to execute this file (chmod +x OMMautomate.sh)."
  echo ""
  echo "Options:"
  echo "  -h, --help: Show this help message and exit."
  echo "  -nt: Number of threads to use in GROMACS simulations (default is 1)."
  echo ""
  echo "Prerequisites:"
  echo "  - .mdp files are located in the home directory."
  echo "  - The PDB files to be used are downloaded in the working directory."
  echo "  - The DSSP library is located in /usr/bin/dssp."
  echo "  - insane executable, Martinate.py, and all required .itp files are in the home directory."
  echo ""
}

# Parse command line arguments
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
  help
  exit 0
fi

# Set default number of threads
nt=1

# Parse optional arguments
while [[ $# -gt 1 ]]
do
key="$2"

case $key in
    -nt|--threads)
    nt="$3"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    echo "Unknown option: $2"
    help
    exit 1
    ;;
esac
done

pdb_code=$1
cg_pdb=${pdb_code}-cg.pdb
cg_top=${pdb_code}-cg.top

# Check if PDB code is provided as an argument
if [ -z "$1" ]
  then
    echo "Please provide a PDB code as an argument"
    exit 1
fi

# Coarse-graining
../martinize.py -f ${pdb_code}.pdb -o ${cg_top} -x ${cg_pdb} -dssp /usr/bin/dssp -p backbone -ff martini22

# Building initial configuration
../insane -u POPC:5.5 -u CHOL:0.5 -u SAPE:4 -alname SAPE -alhead 'E P' -allink 'G G' -altail 'DDDDC CCCC' -l POPC:5.5 -l CHOL:0.5 -l PAPI:2 -l SAPE:2 -alname SAPE -alhead 'E P' -allink 'G G' -altail 'DDDDC CCCC' -d 10 -o system.gro -p topol.top -f ${cg_pdb} -center -pbc hex -sol W -salt 0 -excl -1
