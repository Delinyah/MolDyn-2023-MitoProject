# Lipid fingerprints
# Here, for every frame, the particle positions in the XY plane are binned
# and for each bin the distance to the nearest bin containing protein particles is determined,
# together with the counts of particles of each lipid type per bin.
# From this the cumulative particle counts are determined per frame as a function of distance.
# The raw counts are divided by the numbers of particles per lipid to get
# the cumulative number of lipids of each type as a function of distance.
# These are then normalized to have a sum of one, yielding composition vectors,
# which can be regarded as fingerprints for the region within the specified distance.
# At the maximum distance, the composition is the overall composition
# and this can be used to determine relative compositions allowing assessment of enrichment and depletion.
# The relative compositions are presented on a $\log_{1.1}$ scale,
# such that 1 and -1 indicate 10% enrichment and depletion, respectively.
# This should suffice for the purpose here.
# It is noted that further analysis can take the covariance matrix of the composition vector into account,
# from which it could be seen whether specific lipid types are coupled
# (e.g., negatively charged lipids might exchange and thus be negatively coupled in the covariance matrix).

################

print("Importing modules")

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
import sys
import seaborn as sns

################
gro_file = sys.argv[1]
xtc_file = sys.argv[2]
output_folder = sys.argv[3]

u = mda.Universe(gro_file, xtc_file)
binsize = 2 # Angstrom
nbins = np.round(u.trajectory.ts.dimensions / binsize).astype(int)
binsize = u.trajectory.ts.dimensions / nbins

print("Created your universe and set standard binning settings")

################



# selections
protein = u.select_atoms('protein') # Should select TMD...
membrane = u.select_atoms('not protein')
linkers = membrane.select_atoms('name AM1 or name AM2 or name GL1 or name GL2')

# Function to get the atomgroups for each lipid and each leaflet
#Function to get the upper and lower leaflet
def get_leaflet(universe, headgroup='PO4'):
    #Select the headgroups
    headgroups = universe.select_atoms('name ' + headgroup)
    #Select residues that are lipids, they are not water and they are not ions
    lipids = universe.select_atoms('not protein and not type W and not type NA* and not type CL*')
    #Calculate the center of mass of the headgroups
    center = headgroups.center_of_mass()
    #The lower leaflet is the one that is below the center of mass
    lower_leaflet_headgroups = headgroups.select_atoms('prop z < ' + str(center[2]))
    #The upper leaflet is the one that is above the center of mass
    upper_leaflet_headgroups = headgroups.select_atoms('prop z > ' + str(center[2]))
    #Get the resids of the headgroups in the upper and lower leaflet
    upper_leaflet_resids = upper_leaflet_headgroups.resids
    lower_leaflet_resids = lower_leaflet_headgroups.resids
    #Select the lipids in the upper and lower leaflet
    upper_leaflet = lipids.select_atoms('resid ' + ' '.join(str(x) for x in upper_leaflet_resids))
    lower_leaflet = lipids.select_atoms('resid ' + ' '.join(str(x) for x in lower_leaflet_resids))

    #Make an atomgroup for each unique lipid type and save it in a list
    lipids_upper = []
    lipids_lower = []
    unique_lipid_upper = np.unique(upper_leaflet.resnames)
    unique_lipid_lower = np.unique(lower_leaflet.resnames)
    for lipid in unique_lipid_upper:
        lipids_upper.append(upper_leaflet.select_atoms('resname ' + lipid))
    for lipid in unique_lipid_lower:
        lipids_lower.append(lower_leaflet.select_atoms('resname ' + lipid))
    #Join the lipids in the upper and lower leaflet in a list
    lipids = lipids_upper + lipids_lower

    '''Start: This section can be deleted'''
    #Double check that the number of residues in the upper and lower leaflet is as expected in the .top
    #Get the types of residues and the number of residues in each leaflet
    #Create a dictionary to store the number of residues of each type
    ul_residues = {}
    ll_residues = {}
    #Get the types of residues in the upper leaflet
    for res in upper_leaflet.residues:
        if res.resname in ul_residues:
            ul_residues[res.resname] += 1
        else:
            ul_residues[res.resname] = 1
    #Get the types of residues in the lower leaflet
    for res in lower_leaflet.residues:
        if res.resname in ll_residues:
            ll_residues[res.resname] += 1
        else:
            ll_residues[res.resname] = 1
    #Print the number of residues of each type in each leaflet
    print("The upper and lower leaflet were separated in the following way:")
    print('Number of residues in the upper leaflet:')
    for key in ul_residues:
        print('{}: {}'.format(key, ul_residues[key]))
    print('Number of residues in the lower leaflet:')
    for key in ll_residues:
        print('{}: {}'.format(key, ll_residues[key]))
    '''End: This section can be deleted'''

    return lipids_lower, lipids_upper, lipids

#Get the atomgroups for each lipid and each leaflet
lipids_lower, lipids_upper, lipids = get_leaflet(membrane)
membrane.select_atoms('resname W')

#This line works if we do not want to separete the lipids in the upper and lower leaflet
#lipids = [ u.select_atoms('resname ' + lip) for lip in set(u.select_atoms('not protein').resnames) ]
pangles = protein.positions[:,  2] * (2 * np.pi / membrane.dimensions[2])
mangles = membrane.positions[:,  2] * (2 * np.pi / membrane.dimensions[2])
langles = linkers.positions[:, 2] * (2 * np.pi / membrane.dimensions[2])

linz = np.arctan2(np.sin(langles).mean(), np.cos(langles).mean())

langles += np.pi - linz
langles %= 2 * np.pi
langles -= np.pi

mangles += np.pi - linz
mangles %= 2 * np.pi
mangles -= np.pi

pangles += np.pi - linz
pangles %= 2 * np.pi
pangles -= np.pi

up = langles > 0

uplim = np.quantile(langles[up], 0.9)
lolim = np.quantile(langles[~up], 0.1)

tmd = np.where((pangles > lolim) & (pangles < uplim))[0]

plt.scatter(membrane.positions[:, 0], mangles, s=1, c='#999999', label='Solvent')
plt.scatter(protein.positions[:, 0], pangles, s=1, c='#000099', label='Rest of Protein')
plt.scatter(protein.positions[tmd, 0], pangles[tmd], s=5, c='red', zorder=10, label='TMD Protein')

plt.scatter(linkers.positions[up, 0], langles[up], s=5, c='orange', label='Lipid Linkers Upper Leaflet')
plt.scatter(linkers.positions[~up, 0], langles[~up], s=5, c='orange', label='Lipid Linkers Lower Leaflet')

plt.axhline(0, linewidth=3, c='k')
plt.axhline(uplim, linewidth=2, c='k')
plt.axhline(lolim, linewidth=2, c='k')

plt.xlabel('x-position (nm)')
plt.ylabel('Angle (radians)')
plt.title('Membrane and protein angles')

plt.legend(fontsize='small')

plt.savefig(f"{output_folder}/Membrane_Protein_Angles", dpi=400)

protein = protein.indices[tmd]

print("Done with main selections and plotted the system")

##########################

print("Now defining the process_frame function")

def process_frame(frame, protein, lipids, nbins, binsize, distbins):
    pos = frame.positions
    box = frame.dimensions[:3]

    xybins = ((pos[:, :2] % box[:2]) * (nbins[:2] / box[:2])).astype(int)

    # Cells occupied by protein
    pbins = xybins[protein]
    pcount = np.bincount(pbins[:, 0] * nbins[1] + pbins[:, 1])
    pwhere = np.where(pcount)[0]
    pcells = np.stack((pwhere // nbins[1], pwhere % nbins[1]), axis=1)
    # Distances of cells from protein
    cells = np.mgrid[:nbins[0], :nbins[1]].T.reshape((-1, 2))
    distances = (((cells[:, None, :] - pcells[None, :]) * binsize)**2).sum(axis=2).min(axis=1) ** 0.5
    
    # Bincounts per lipid type
    counts = {}
    for lip in lipids:
        lbins = xybins[lip.indices]
        lcount = np.bincount(lbins[:, 0] * nbins[1] + lbins[:, 1], minlength=cells.shape[0])
        counts[lip] = lcount
        
    # Composition as function of distance
    lo = -1
    fpr = []
    for hi in distbins:
        m = (distances > lo) & (distances <= hi)
        fpr.append([c[m].sum() for l, c in counts.items()])
        lo = hi
    return fpr
print("Defined process_frame function")

##########################

print("Now I will process your frames. This may take a while.")
# For 5000 frames this takes a few minutes on a MacBook Pro.
# per frame per distance a cumulative number composition vector
# frames * distancebins * lipidtypes
start_time = time.time() # record start time
fingerprints = []
fingerprints = []
distbins = np.linspace(0, 80, 81)
for frame in u.trajectory[::10]:
    fp = process_frame(frame, protein, lipids, nbins[:2], binsize[:2], distbins)
    fingerprints.append(fp)
end_time = time.time() # record end time
elapsed = end_time - start_time
print(f"Processed frames in {elapsed:.2f} seconds")
names_upper = [lip[0].resname + '_u' for lip in lipids_upper]
names_lower = [lip[0].resname + '_l' for lip in lipids_lower]
names = names_upper + names_lower
# Indices for upper and lower leaflets
upper_indices = [i for i, name in enumerate(names) if name.endswith('_u')]
lower_indices = [i for i, name in enumerate(names) if name.endswith('_l')]
#This line works if we do not want to separete the lipids in the upper and lower leaflet
#names = [lip[0].resname for lip in lipids]
sizes = [ len(l.split('residue')[0]) for l in lipids ]
F = np.array(fingerprints).cumsum(axis=1)
M = F.mean(axis=0) / sizes
P = M / M.sum(axis=1)[:, None]
L = np.log(P / P[-1]) / np.log(1.1)
print("Log1.1 enrichment with respect to total as function of distance from protein.")
print("   ", *names, sep="  ")
for d, f in zip(distbins, np.round(L, 2)):
    print(d, f)

print('---')

print("Fingerprint percentages as function of distance from protein.")
print("   ", *names, sep="  ")
for d, f in zip(distbins, np.round(100*P, 2)):
    print(d, f)

##########################

# Plot for upper leaflet
plt.figure()

for i in upper_indices:
    plt.plot(distbins, L[:, i])

plt.legend(names_upper)
plt.xlabel('Distance from protein (...)')
plt.ylabel('log1.1 enrichment')
plt.title('Log1.1 enrichment per lipid type as function of distance from protein (Upper)')
plt.axhline(0,color='k', linestyle='--', linewidth=1)
plt.savefig(f"{output_folder}/Log1-1_enrichment_upper", dpi=400)

# Plot for lower leaflet
plt.figure()

for i in lower_indices:
    plt.plot(distbins, L[:, i])

plt.legend(names_lower)
plt.xlabel('Distance from protein (...)')
plt.ylabel('log1.1 enrichment')
plt.title('Log1.1 enrichment per lipid type as function of distance from protein (Lower)')
plt.axhline(0,color='k', linestyle='--', linewidth=1)
plt.savefig(f"{output_folder}/Log1-1_enrichment_lower", dpi=400)

print("Plotted log1.1 enrichments.")

##########################

composition_matrix_upper = P[:, upper_indices]
composition_matrix_lower = P[:, lower_indices]

cormat_upper = np.corrcoef(composition_matrix_upper.T)
cormat_lower = np.corrcoef(composition_matrix_lower.T)

print("Correlation matrix upper:")
print(cormat_upper)

print("Correlation matrix lower:")
print(cormat_lower)

# Heatmap for upper leaflet
plt.figure(figsize=(4, 4))
sns.heatmap(cormat_upper, annot=False, cmap="coolwarm", xticklabels=names_upper, yticklabels=names_upper)
plt.title("Correlation Matrix of Lipid Compositions (Upper)")
plt.savefig(f"{output_folder}/Cormat_lip_comp_upper", dpi=400)

# Heatmap for lower leaflet
plt.figure(figsize=(4, 4))
sns.heatmap(cormat_lower, annot=False, cmap="coolwarm", xticklabels=names_lower, yticklabels=names_lower)
plt.title("Correlation Matrix of Lipid Compositions (Lower)")
plt.savefig(f"{output_folder}/Cormat_lip_comp_lower", dpi=400)

print("done with this script")
