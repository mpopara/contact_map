# contact_map
Python script for mapping the presence of inter-residue contacts in a biomolecule

## General description

This script analyses the ensemble of biomolecular structures with corresponding population fractions (weights), in terms of intramolecular interactions.
For a selected type of atoms, e.g. backbone N atoms, calculated are distances for each conformer in the ensemble and for each pair of atoms.
Choosing backbone N atoms for this analysis is particularly convenient because it enables comparison with NMR data obtained for <sup>15</sup>N-labelled protein.
Chemical Shift Perturbations (CSPs) obtained in NMR experiment are indicative of changes in local magnetic field of NMR active atoms (here <sup>15</sup>N),
which in turn reflect occurence of inter-residue contacts. Once CSPs are mapped on a protein structure, areas which experience stronger perturbations
can be easily identified and compared with computational models, obtained from e.g. all-atom MD simulations.

Based on a user-defined cut-off value, distances between all pair-wise combinations of N atoms are converted to binary information: 1-contact occurs, 0-no contact. 
Choice of cut-off value depends on the type of the experimental data to which results of this analysis are to be compared, and type of atoms for which distances are calculated.
In the next step, using the weights of the ensemble members, performed is the weighted averaging of the contact maps.
Lastly, ensemble-averaged contact map is visualized and saved as png and svg file. 

![contact_map](https://github.com/mpopara/contact_map/assets/40856779/91d919ae-5d2e-42c3-91fa-9c84da363ca4)


## Example data and input file requirements

In the folder example_data, provided are exemplary input files: 
* ensemble of conformational models as a single trajectory file in .dcd format (any mdtraj-compatible trajectory format 
will work as well- .xtc, .nc, .h5, .pdb...). While the origin of the conformational ensemble and population fractions is irrelevant, it is necessary that conformational models are all-atom, 
since the goal of the script is to identify presence of contacts on residue/atomistic level.
* topology file in .pdb format
* .dat file contaning weights (population fractions) of ensemble members. This file is of a size N<sub>conformers</sub> x 2, where the first column contains indices of the ensemble members, and second columns their corresponding weights.
  Order or ensemble members in trajectory file should be the same as in the weights file. 


## System requirements

contacts_map.py is a python script built on Python 3.8.8. Script was tested under the following configuration:

* Windows 10
* Python 3.8.8
* numpy 1.23.0
* mdtraj 1.9.4
* matplotlib 3.7.1

## Authors

* Milana Popara
