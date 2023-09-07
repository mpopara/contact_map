import numpy as np
import mdtraj as md
import itertools
from itertools import *
import matplotlib.pyplot as plt

##### input parameters and data paths #######

sel_1 = "name N"
sel_2 = "name N"

threshold = 12 # distance threshold in Angstrom


traj_path = 'C:/user/folder/trajectory.dcd' # path to the trajectory file
top_path = 'C:/user/folder/topology.pdb' # path to the topology file

weights_path = 'C:/user/folder/conformer_weights.dat' # path to the weights file of ensemble members


save_path = 'C:/user/folder/' # path where results will be saved


######## loading input files ########################

w = np.loadtxt(weights_path, delimiter=' ',usecols=[1])
w=np.atleast_2d(w).T
	
traj=md.load(traj_path,top=top_path)
traj.unitcell_angles = None
traj.unitcell_lengths = None
traj.unitcell_vectors = None	


assert len(w)==len(traj), "Size of weight file does not match number of conformers in trajectory" # make sure than n_conformers in weights file matches n_conformers from trajectory


###################### compute distance map ###########################

idx_1 = traj.topology.select(sel_1)
idx_2 = traj.topology.select(sel_2)
idx_pairs = list(itertools.product(idx_1, idx_2)) # makes list of all pair-wise combinations of idx_1 and idx_2 indices: 0-0,0-1,0-2,..,1-0,1-1,1-2,..

distances = md.compute_distances(traj, idx_pairs)*10 # native mdtraj unit is nm

cm=(distances<threshold).astype(int) # convert distances to binary information/contact map

weighted_cm=np.sum(w*cm,axis=0)/np.sum(w)*100 #ensemble weighted average of contact maps


#################### visualization ############################

resid_1=np.arange(1,traj.n_residues+1)
resid_2=np.arange(1,traj.n_residues+1)
resid_pairs = list(itertools.product(resid_1, resid_2)) #for visualization, create pairs of residue numbers,
                                                        # instead of pairs of N indices

x_val=[x[0] for x in resid_pairs]
y_val=[x[1] for x in resid_pairs]


fig, ax = plt.subplots()
im=ax.scatter(x_val, y_val, s=5, c=weighted_cm, cmap='viridis')
cbar=plt.colorbar(im).set_label(label='contact occurence [%]', size=12)

ax.set_xlabel('residue number')
ax.set_ylabel('residue number')

ax.set_xlim(1,traj.n_residues+1)
ax.set_ylim(1,traj.n_residues+1)

plt.tight_layout()
plt.savefig(save_path+'contact_map.png', dpi=300, transparent=True)
plt.savefig(save_path+'contact_map.svg', transparent=True)
plt.close(fig)


