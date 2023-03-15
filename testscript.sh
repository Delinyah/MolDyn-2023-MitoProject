#!/bin/bash

# Help function 
function help {
  echo "" 
  echo -e "\033[38;5;226mAutomated Workflow for Outer Mitochondrial Membrane Simulation Modeling\033[0m"   
  echo -e "\033[38;5;208mTitle: mainone.sh (for 'MArtinize-INsanify Outer-membraNE...'\033[0m"
  echo -e "\033[38;5;34mAuthor: Delinyah C. Koning (as of March 9, 2023)\033[0m"
  echo -e "\033[38;5;34mUniversity of Groningen, FSE Faculty, GBB Institute, Molecular Dynamics Group (2023)\033[0m" 
  echo ""
  echo "Usage: ./mainone.sh <pdb_code> [-nt <number_of_threads>]"
  echo "Don't forget to give permission to execute this file (chmod +x mainone.sh)."
  echo ""
  echo -e "\033[38;5;226mWithout editing the parameter files, the run is executed as follows:\033[0m"
  echo "Martinize coarse-grains protein"
  echo "Insane builds a charge-neutralized (NaCl) coarse-grained system with -excl set to -1 so that water can be placed everywhere (e.g. also inside barrel proteins)."
  echo "1000-step EM"
  echo "1500ps NVT (dt 0.03, nstxout 20, v-rescale coupling system)"
  echo "1500ps NPT (dt 0.03, nstxout 20, parrinello-rahman pressure coupling, semiisotropic)"
  echo "30ns production run (dt 0.03, nstep 1000000, nstxout 200)" 
  echo ""
  echo "Notes:"  
  echo "The pdb-code must be specified without .pdb extension (just the 4-letter code)."
  echo "The pdb-file is the only file that should be present in the working directory. Fetching structures from PDB or OPM will be available in the future."
  echo "After running, all files in the working directory will be moved (not copied) to a folder that is named after the pdb-code."
  echo "Moving files will be done by moving all (*) files with a certain extension type. The .pdb files are specified manually to prevent mixup."
  echo ""
  echo "Flag options:"
  echo "  -h, --help: Show this help message and exit."
  echo "  -nt: Number of threads to use in GROMACS simulations (default is 1)."
  echo ""
  echo "Prerequisites:"
  echo "  - .mdp files are located in the home directory (they will be sourced and copied to the working directory)."
  echo "  - The PDB files to be used MUST be downloaded in the working directory."
  echo "  - The DSSP library is located in /usr/bin/dssp. If not, adapt script."
  echo "  - Insane executable, Martinate.py, and all required .itp files are in the home directory. They will not be copied to the working directory."
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

#Cleanup
echo -e "\033[38;5;34mCreating an output folder...\033[0m"
mkdir ./"${pdb_code}"
mv ${pdb_code}.pdb ./"$1"
cd ./"$1"

# Check if PDB code is provided as an argument
if [ -z "$1" ]
  then
    echo "Please provide a PDB code as an argument"
    exit 1
fi

# Coarse-graining
echo -e "\033[38;5;226mCoarse graining your system...\033[0m"
~/Project/OMM/martinize.py -f ${pdb_code}.pdb -o topol.top -x ${cg_pdb} -dssp /usr/bin/dssp -p backbone -ff martini22

# Building initial configuration
echo -e "\033[38;5;226mHOld on, building system...\033[0m"
~/Project/OMM/insane -u POPC:5.5 -u CHOL:0.5 -u SAPE:4 -alname SAPE -alhead 'E P' -allink 'G G' -altail 'DDDDC CCCC' -l POPC:5.5 -l CHOL:0.5 -l PAPI:2 -l SAPE:2 -alname SAPE -alhead 'E P' -allink 'G G' -altail 'DDDDC CCCC' -d 10 -o system.gro -p topol.top -f ${cg_pdb} -center -pbc hex -sol W -salt 0 -excl -1

# Modify topol.top include statements
echo ';' >> topol.top
~/Project/OMM/insane -u POPC:5.5 -u CHOL:0.5 -u SAPE:4 -alname SAPE -alhead 'E P' -allink 'G G' -altail 'DDDDC CCCC' -l POPC:5.5 -l CHOL:0.5 -l PAPI:2 -l SAPE:2 -alname SAPE -alhead 'E P' -allink 'G G' -altail 'DDDDC CCCC' -d 10 -o system.gro -p topol.top -f ${cg_pdb} -center -pbc hex -sol W -salt 0 -excl -1 2>&1 | tee -a topol.top

# Add needed include topology statements to topol.top
sed -i 's/#include "martini.itp"/#include "..\/martini_v2.2.itp"\n#include "..\/SAPE.itp"\n#include "..\/martini_v2.0_ions.itp"\n#include "..\/martini_v2.0_lipids_all_201506.itp"/; s/\bProtein\b/Protein/g' topol.top

# EM
echo -e "\033[38;5;226mEnergy minimizing...\033[0m"
echo "${min_mdp}" > minimization.mdp
gmx grompp -p topol.top -f ~/Project/OMM/minimization.mdp -c system.gro -o minimization.tpr -maxwarn 1
gmx mdrun -v -deffnm em -s minimization.tpr -nt $nt

# NVT
echo -e "\033[38;5;226mNVT equilibration...\033[0m"
echo "${nvt_mdp}" > nvt.mdp
gmx grompp -f ~/Project/OMM/nvt.mdp -c em.gro -p topol.top -o nvt.tpr -maxwarn 1
gmx mdrun -v -deffnm nvt -s nvt.tpr -nt $nt

# NPT
echo -e "\033[38;5;226mNPT equilibration...\033[0m"
echo "${npt_mdp}" > npt.mdp
gmx grompp -f ~/Project/OMM/npt.mdp -c nvt.gro -p topol.top -o npt.tpr -maxwarn 1
gmx mdrun -v -deffnm npt -s npt.tpr -nt $nt

# Production run
echo -e "\033[38;5;226mNow running the real deal...\033[0m"
echo "${run_mdp}" > run.mdp
gmx grompp -f ~/Project/OMM/run.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 1
gmx mdrun -v -deffnm md -s md.tpr -nt $nt

#Analysis (echoes '0' for whole system)
echo -e "\033[38;5;226mCalculating RMSD, RMSF, and Rg after run...\033[0m"
yes 0 | head -n 2 | gmx rms -s md.tpr -f md.xtc -o rmsd_md.xvg
echo 0 | gmx rmsf -s md.tpr -f md.xtc -o rmsf_md.xvg
echo 0 | gmx gyrate -s md.tpr -f md.xtc -o gyrate_md.xvg

#Shoutouts
echo " "
echo -e "\033[38;5;208m'The computer was born to solve problems that did not exist before.' — Bill Gates, Microsoft founder and former CEO\033[0m"
echo " "
echo -e "\033[38;5;226mYou made it to the end of the script. That can't be right...\033[0m"
echo -e "\033[38;5;226m:)\033[0m"
echo " "
