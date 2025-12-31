import numpy as np
from itertools import product

def fetch_supercell_motif(motif: np.ndarray, supercell_dims: np.ndarray):
    
    Nx,Ny,Nz = supercell_dims
    
    supercell_motif = []
    for i,j,k in product(range(0,Nx),range(0,Ny),range(0,Nz)):
        for pos in motif:
            new_pos = pos + np.array([i,j,k])
            supercell_motif.append(new_pos/supercell_dims)
    
    return np.array(supercell_motif)
            

if __name__ == '__main__':

    motif = np.array([
 [0.000000 ,  0.000000,   0.000000],
 [0.500000 ,  0.500000,   0.000000],
 [0.500000 ,  0.000000 ,  0.500000],
 [0.000000 ,  0.500000 ,  0.500000],
 [0.250000 ,  0.250000 ,  0.250000],
 [0.750000 ,  0.250000 ,  0.250000],
 [0.250000 ,  0.750000 ,  0.250000],
 [0.750000 ,  0.750000 ,  0.250000],
 [0.250000  , 0.250000 ,  0.750000],
 [0.750000  , 0.250000 ,  0.750000],
 [0.250000 ,  0.750000 ,  0.750000],
 [0.750000 ,  0.750000 ,  0.750000],
 ])
    supercell_dims = np.array([2,1,1])
    
    print(fetch_supercell_motif(motif,supercell_dims))
    


            