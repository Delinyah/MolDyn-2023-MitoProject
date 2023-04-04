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

################
gro_file = sys.argv[1]
xtc_file = sys.argv[2]
output_folder = sys.argv[3]

u = mda.Universe(gro_file, xtc_file)
binsize = 2 # Angstrom
nbins = np.round(u.trajectory.ts.dimensions / binsize).astype(int)
binsize = u.trajectory.ts.dimensions / nbins
################

print("Created your universe and set standard binning settings")

# selections
protein = u.select_atoms('protein') # Should select TMD...
membrane = u.select_atoms('not protein')
linkers = membrane.select_atoms('name AM1 or name AM2 or name GL1 or name GL2')

# Function to get the atomgroups for each lipid and each leaflet
def get_leaflet(universe, headgroup='PO4'):
    
    headgroups = universe.select_atoms('name ' + headgroup)
    lipids = universe.select_atoms('not protein and not type W and not type NA* and not type CL*')
    center = headgroups.center_of_mass()
    
    lower_leaflet_headgroups = headgroups.select_atoms(f'prop z < {center[2]}')
    upper_leaflet_headgroups = headgroups.select_atoms(f'prop z > {center[2]}')
    
    upper_leaflet_resids = upper_leaflet_headgroups.resids
    lower_leaflet_resids = lower_leaflet_headgroups.resids
    
    upper_leaflet = lipids.select_atoms('resid ' + ' '.join(str(x) for x in upper_leaflet_resids))
    lower_leaflet = lipids.select_atoms('resid ' + ' '.join(str(x) for x in lower_leaflet_resids))
    
    unique_lipid_upper = np.unique(upper_leaflet.resnames)
    unique_lipid_lower = np.unique(lower_leaflet.resnames)
    lipids_upper = [upper_leaflet.select_atoms(f'resname {lipid}') for lipid in unique_lipid_upper]
    lipids_lower = [lower_leaflet.select_atoms(f'resname {lipid}') for lipid in unique_lipid_lower]
    
    lipids = lipids_upper + lipids_lower
    
    return lipids_lower, lipids_upper, lipids

#Get the atomgroups for each lipid and each leaflet
lipids_lower, lipids_upper, lipids = get_leaflet(membrane)

#This hashed line below works if we do not want to separate the lipids in the upper and lower leaflet
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

plt.savefig("Membrane_Protein_Angles", dpi=400)

plt.show()

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
fingerprints = []

distbins = np.linspace(0, 80, 81)

start_time = time.time() # record start time

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

plt.figure()

for i in upper_indices:
    plt.plot(distbins, L[:, i])

plt.legend(names_upper)
plt.xlabel('Distance from protein (...)')
plt.ylabel('log1.1 enrichment')
plt.title('Log1.1 lipid enrichment in upper leaflet as function of distance from protein')
plt.axhline(0,color='k', linestyle='--', linewidth=1)
plt.savefig("Log1-1_enrichment_upper", dpi=400)

# Plot for lower leaflet
plt.figure()

for i in lower_indices:
    plt.plot(distbins, L[:, i])

plt.legend(names_lower)
plt.xlabel('Distance from protein (...)')
plt.ylabel('log1.1 enrichment')
plt.title('Log1.1 lipid enrichment in lower leaflet as function of distance from protein')
plt.axhline(0,color='k', linestyle='--', linewidth=1)
plt.savefig("Log1-1_enrichment_lower", dpi=400)

print("Plotted log1.1 enrichments.")

##########################

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Visualization
# Slice the data based on upper and lower leaflet indices
L_upper = L[:, upper_indices]
L_lower = L[:, lower_indices]

# Heatmap for upper leaflet
fig, ax = plt.subplots(figsize=(8, 4))
sns.heatmap(L_upper.T, cmap='coolwarm', cbar_kws={'label': 'Log1.1 Enrichment'}, xticklabels=5, yticklabels=names_upper, ax=ax)
ax.set_xlabel('Distance from Protein')
ax.set_ylabel('Lipid Types')
ax.set_title('Log1.1 Enrichment with respect to Total as Function of Distance from Protein; upper leaflet')
plt.savefig("Log1-1_enrichment_upper_hm", dpi=400)

# Heatmap for lower leaflet
fig, ax = plt.subplots(figsize=(8, 4))
sns.heatmap(L_lower.T, cmap='coolwarm', cbar_kws={'label': 'Log1.1 Enrichment'}, xticklabels=5, yticklabels=names_lower, ax=ax)
ax.set_xlabel('Distance from Protein')
ax.set_ylabel('Lipid Types')
ax.set_title('Log1.1 Enrichment with respect to Total as Function of Distance from Protein; lower leaflet')
plt.savefig("Log1-1_enrichment_lower_hm", dpi=400)

print("Plotted log1.1 enrichments as heatmap")

##########################

#correlation matrix
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
plt.savefig("Cormat_lip_comp_upper", dpi=400)

# Heatmap for lower leaflet
plt.figure(figsize=(4, 4))
sns.heatmap(cormat_lower, annot=False, cmap="coolwarm", xticklabels=names_lower, yticklabels=names_lower)
plt.title("Correlation Matrix of Lipid Compositions (Lower)")
plt.savefig("Cormat_lip_comp_lower", dpi=400)

print("Plotted correlation matrix")
print("Done with this .py script")
