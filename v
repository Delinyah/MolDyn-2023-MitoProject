#Analysis
gmx rms -s npt.tpr -f npt.xtc -o rmsd_npt.xvg
gmx rmsf -s npt.tpr -f npt.xtc -o rmsf_npt.xvg
gmx gyrate -s npt.tpr -f npt.xtc -o gyrate_npt.xvg
