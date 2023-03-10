#make selections
select POPC, resn POPC
select CHOL, resn CHOL
select POPE, resn POPE
select SAPE, resn SAPE
select PAPI, resn PAPI
select CDL1, resn CDL1
select water, resn H
select protein, resn ala+arg+asn+asp+cys+gln+glu+gly+his+ile+leu+lys+met+phe+pro+ser+thr+trp+tyr+val

#color by element
color slate, (name C* and POPC)
color cyan, (name C* and CHOL)
color green, (name C* and POPE)
color orange, (name C* and SAPE)
color magenta, (name C* and PAPI)
color yellow, (name C* and CDL1)

#load trajectory (assumes file output is prefixed by 'npt' during mdrun; if not, change script)
load npt.xvg, npt
