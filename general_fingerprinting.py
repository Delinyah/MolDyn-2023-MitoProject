import argparse
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time

parser = argparse.ArgumentParser(description='Fingerprinting script for MD trajectories')
parser.add_argument('gro_file', type=str, help='Input .gro file')
parser.add_argument('xtc_file', type=str, help='Input .xtc file')

args = parser.parse_args()

gro_file = args.gro_file
xtc_file = args.xtc_file

binsize = 2  # Angstrom

# Replace the hard-coded file paths with the gro_file and xtc_file variables
u = mda.Universe(gro_file, xtc_file)
nbins = np.round(u.trajectory.ts.dimensions / binsize).astype(int)
binsize = u.trajectory.ts.dimensions / nbins

## selections

protein = u.select_atoms('protein')  # Should select TMD...
membrane = u.select_atoms('not protein')
linkers = membrane.select_atoms('name AM1 or name AM2 or name GL1 or name GL2')
lipids = [u.select_atoms('resname ' + lip) for lip in set(u.select_atoms('not protein').resnames)]

pangles = protein.positions[:, 2] * (2 * np.pi / membrane.dimensions[2])
mangles = membrane.positions[:, 2] * (2 * np.pi / membrane.dimensions[2])
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

plt.show()

protein = protein.indices[tmd]

## frame processing

def process_frame(frame, protein, lipids, nbins, binsize, distbins):
    pos = frame.positions
    box = frame.dimensions[:3]
    box_half = box / 2

    dist = np.zeros((len(protein), len(lipids)))
    for i, p in enumerate(protein):
        for j, l in enumerate(lipids):
            lpos = pos[l.indices]
            diff = lpos - pos[p]
            diff -= box * np.round(diff / box)
            dist[i, j] = np.sqrt((diff * diff).sum(axis=1)).min()

    hist, _ = np.histogram(dist, bins=distbins)
    return hist / len(protein)

## main loop

start_time = time.time()
hist = np.zeros((len(lipids), nbins[0]))

distbins = np.linspace(0, nbins[0] * binsize[0], nbins[0] + 1)

for frame in u.trajectory:
    hist += process_frame(frame, protein, lipids, nbins, binsize, distbins)
    print(f"Processed frame {frame.frame} / {len(u.trajectory)}")

hist /= len(u.trajectory)

## output

df = pd.DataFrame(hist, columns=distbins[:-1], index=[lip.resname for lip in lipids])
df.to_csv('lipid_fingerprint.csv')

print(f"Done! Time taken: {time.time() - start_time:.2f} seconds")
