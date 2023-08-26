#!/bin/bash

# Help function 
function help {
  echo "" 
  echo -e "\033[38;5;226mAutomated Workflow for Outer Mitochondrial Membrane Simulation Modeling\033[0m"   
  echo -e "\033[38;5;208mTitle: mainone.sh (for 'MArtinize-INsanify Outer-membraNE...'\033[0m"
  echo -e "\033[38;5;34mAuthor: Delinyah C. Koning (last edited on March 16th, 2023)\033[0m"
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
  echo "by default a 100ns production run (dt 0.03, nstep 3333333, nstxout 1000 >>> 3333 frames total, 30ps per frame)" 
  echo ""
  echo "Notes:"  
  echo "The pdb-code must be specified without .pdb extension (just the 4-letter code)."
  echo "The pdb-file is the only file that should be present in the working directory (+ this .sh script). Fetching structures from PDB or OPM will be available in the future."
  echo "During running, the pdb-input in the initial working directory will be moved (not copied) to a folder that is named after the pdb-code."
  echo "During the script, the new directory will become the new working directory; this makes sure that all output is collected."
  echo "The topol.top file after coarse-graining is adapted by appending it with the stderr output after execution of the insane command. Do NOT change."
  echo ""
  echo "Flag options:"
  echo "  -h, --help: Show this help message and exit."
  echo "  -nt: Number of threads to use in GROMACS simulations (default is 1)."
  echo ""
  echo "Prerequisites:"
  echo "  - .mdp files are located in the home directory (specify path as needed)."
  echo "  - The PDB files to be used MUST be downloaded in the working directory."
  echo "  - The DSSP library is assumed to be located in /usr/bin/dssp. If not, adapt script."
  echo "  - Insane executable, Martinate.py, and all required .itp files are in the home directory. They will not be copied to the working directory."
  echo ""
}

# Error function
function error_exit {
  echo -e "\033[38;5;196mMwah-Mwah..... Something is wrong here...\033[0m" >&2
  exit 1
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
/martini/delinyah/Project/OMM/martinize.py -f ${pdb_code}.pdb -o topol.top -x ${cg_pdb} -dssp /usr/bin/dssp -p backbone -ff martini22 || error_exit

# Make sure that when the following command edits topol.top no new stuff is added to existing lines
echo ';' >> topol.top

# Building initial configuration + writing text file with stderr output
echo -e "\033[38;5;226mHold on, building system...\033[0m"
/martini/delinyah/Project/OMM/insane -u POPC:5.4 -u SAPE:2.0 -alname SAPE -alhead 'E P' -allink 'G G' -altail 'DDDDC CCCC' -u PAPI:2.2 -u POPS:0.3 -u POPA:0.1 -l POPC:5.4 -l SAPE:3.8 -alname SAPE -alhead 'E P' -allink 'G G' -altail 'DDDDC CCCC' -l PAPI:0.4 -l POPS:0.1 -l POPA:0.1 -alname POPA -alhead 'P' -allink 'G G' -altail 'CDCC CCCC' -l CDL1:0.2 -d 10 -o system.gro -f ${cg_pdb} -center -orient -pbc hex -sol W -salt 0 -excl -1 2>&1 | tee -a topol.top || error_exit

# Add other needed include topology statements to topol.top
sed -i 's/#include "martini.itp"/#include "..\/martini_v2.2.itp"\n#include "..\/SAPE.itp"\n#include "..\/martini_v2.0_ions.itp"\n#include "..\/martini_v2.0_lipids_all_201506.itp"/; s/\bProtein\b/Protein/g' topol.top

# EM
echo -e "\033[38;5;226mEnergy minimizing...\033[0m"
echo "${min_mdp}" > minimization.mdp
gmx grompp -p topol.top -f /martini/delinyah/Project/OMM/minimization.mdp -c system.gro -o minimization.tpr -maxwarn 1 || error_exit
gmx mdrun -v -deffnm em -s minimization.tpr -nt $nt || error_exit

# NVT
echo -e "\033[38;5;226mNVT equilibration...\033[0m"
echo "${nvt_mdp}" > nvt.mdp
gmx grompp -f /martini/delinyah/Project/OMM/nvt.mdp -c em.gro -p topol.top -o nvt.tpr -maxwarn 1 || error_exit
gmx mdrun -v -deffnm nvt -s nvt.tpr -nt $nt || error_exit

# NPT
echo -e "\033[38;5;226mNPT equilibration...\033[0m"
echo "${npt_mdp}" > npt.mdp
gmx grompp -f /martini/delinyah/Project/OMM/npt.mdp -c nvt.gro -p topol.top -o npt.tpr -maxwarn 1 || error_exit
gmx mdrun -v -deffnm npt -s npt.tpr -nt $nt || error_exit

# Production run
echo -e "\033[38;5;226mStarting production run for ${pdb_code}...\033[0m" 
echo "${run_mdp}" > run.mdp
gmx grompp -f /martini/delinyah/Project/OMM/run.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 1 || error_exit
gmx mdrun -v -deffnm md -s md.tpr -nt $nt || error_exit

#Analysis (echoes '0' for whole system)
#echo -e "\033[38;5;226mCalculating RMSD, RMSF, and Rg after run...\033[0m"
#yes 1 | head -n 2 | gmx rms -s "$folder/md.tpr" -f "$folder/md.xtc" -o "$output_folder/rmsd_md_protein.xvg" || error_exit
#echo -e "7\n0" | gmx energy -f "$folder/em.edr" -o "$output_folder/energy_em.xvg"  || error_exit
#echo -e "12\n0" | gmx energy -f "$folder/npt.edr" -o "$output_folder/pressure_npt.xvg" || error_exit
#echo -e "11\n0" | gmx energy -f "$folder/nvt.edr" -o "$output_folder/temperature_nvt.xvg" || error_exit

#Shoutouts
echo " "
echo -e "\033[38;5;208m'The computer was born to solve problems that did not exist before.' — Bill Gates, Microsoft founder and former CEO\033[0m"
echo " "
echo -e "\033[38;5;226mEnd of script reached. Are you still awake? Please collect your complimentary rubber duck for debugging at the exit.\033[0m"
echo -e "\033[38;5;226m:)\033[0m"
echo " "

echo -e "\033[38;5;208mFinished production run for ${pdb_code}...\033[0m"
