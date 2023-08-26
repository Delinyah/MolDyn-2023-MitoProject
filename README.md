# MolDyn 2023-MitoProject
**Molecular Dynamics of Mitochondrial Membrane Proteins in Their Native Lipid Environment**
This repository covers all files that were made and used for high-throughput computational modeling of mitochondrial membrane-protein systems. This project at Molecular Dynamics, University of Groningen, was adequate for fulfillment of the bachelor Biology and Medical Laboratory Research at Hanze University of Applied Sciences.<br>
_Group: Molecular Dynamics, University of Groningen (prof. S.J. Marrink)_<br>
_Supervisors: Dr. Tsjerk A. Wassenaar and Msc. Rubi Garcia-Zarmiento_<br>

Mitochondria are unique cell organelles with distinctive dual lipid bilayer configurations. Especially the complex structure of the inner mitochondrial membrane ties membrane organization together with functionality, but also makes the membrane susceptible to structural instability. Stability is largely faxcilitated by protein-lipid interactions that modulate a synchronized molecular interplay between lateral lipid diffusion, membrane curvature, and protein-membrane assembly. Unfortunately, interactions within distinct membrane organizations are difficult to characterize experimentally. In this study, we introduce a systematic high-throughput pipeline that conducts a bulk series of coarse-grained molecular dynamics Gromacs simulations with the Martini2 force field to characterize proteomic and lipidomic interactions in arbitrarily assembled protein-membrane systems of the human mitochondrion, in which we present methods for initial system configuration, simulation, lipid enrichment analysis and trajectory profile clustering. This resulted in bulk-simulation of 147 unique membrane-protein systems, yielding 1830 stable trajectories that were analyzed through spatial enrichment and hierarchical clustering, revealing 28 distinct membrane- and leaflet-dependent clusters. These clusters have implications for membrane morphology and function and suggest (co-)localization of proteins to similar membrane microenvironments. Collectively, this study provides a framework for the characterization and clustering of mitochondrial membrane proteins based on lipid interactions, where we offer a comprehensive data set for future model upscaling and propose a novel theoretical and methodological pipeline for early-phase molecular dynamics modeling studies of large, complex biomolecular systems.

# Aims of the study
 1. Setup of initial system configurations of a comprehensive selection of mitochondrial membrane proteins
 2. Conducting a systematic high-throughput simulation series while defining correct settings and parameters
 3. Characterization of a bulk trajectory series with a focus on leaflet-specific lipid enrichment/depletion profiles
 4. Developing and conducting a trajectory profile clustering method
 5. Proposing a novel pipeline for developing new, detailed computational models of large biomolecular and structurally complex systems

# Contents of this repository
 - Protein selection based on MitoCarta3.0 database by Broad Institute of MIT and Harvard
 - Coarse-grained lipid structures and topologies used for custom membrane assembly
 - Simulation parameters
 - Automation and parallelization scripts
 - Python script for lipid enrichment/depletion analysis (enrichmentsmitochondrion.py)
 - Jupyter Notebook with example for enrichment/depletion analysis (enrichmentsmitochondrion.ipynb)
 - Jupyter Notebook as conducted with example for trajectory profile clustering (bulktrajectorygrouping.ipynb)
