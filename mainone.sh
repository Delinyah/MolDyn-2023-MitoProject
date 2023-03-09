#!/bin/bash

# Help function 
function help {
  echo "" 
  echo -e "\033[38;5;208m'Automated Workflow for Outer Mitochondrial Membrane Simulation Modeling\033[0m"   
  echo -e "\033[38;5;226mTitle: mainone (for MArtinize-INsanify Outer-membraNE...\033[0m"  
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

# Modify topol.top include statements
sed -i 's/#include "martini.itp"/#include "..\/martini_v2.2.itp"\n#include "..\/SAPE.itp"\n#include "Protein_A.itp"\n#include "..\/martini_v2.0_ions.itp"\n#include "..\/martini_v2.0_lipids_all_201506.itp"/; s/\bProtein\b/Protein_A/g' topol.top

# EM
echo "${min_mdp}" > minimization.mdp
gmx grompp -p topol.top -f ../minimization.mdp -c system.gro -o minimization.tpr -maxwarn 1
gmx mdrun -v -deffnm em -s minimization.tpr -nt $nt

# NVT
echo "${nvt_mdp}" > nvt.mdp
gmx grompp -f ../nvt.mdp -c em.gro -p topol.top -o nvt.tpr -maxwarn 1
echo " "gmx mdrun -v -deffnm nvt -s nvt.tpr -nt $nt

# NPT
echo "${npt_mdp}" > npt.mdp
gmx grompp -f ../npt.mdp -c em.gro -p topol.top -o npt.tpr -maxwarn 1
gmx mdrun -v -deffnm npt -s npt.tpr -nt $nt

#Analysis (echoes '0' for whole system)
yes 0 | head -n 2 | gmx rms -s npt.tpr -f npt.xtc -o rmsd_npt.xvg
echo 0 | gmx rmsf -s npt.tpr -f npt.xtc -o rmsf_npt.xvg
echo 0 | gmx gyrate -s npt.tpr -f npt.xtc -o gyrate_npt.xvg

#Shoutouts
echo " "
echo -e "\033[38;5;208m'The computer was born to solve problems that did not exist before.' — Bill Gates, Microsoft founder and former CEO\033[0m"
echo " "
echo -e "\033[38;5;226mYou made it to the end of the script. That can't be right...\033[0m"
echo " "
echo " "
echo " "
