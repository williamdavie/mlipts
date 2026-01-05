import numpy as np
from itertools import product
from ase import Atoms

def fetch_supercell_motif(motif: np.ndarray, supercell_dims: np.ndarray):
    
    Nx,Ny,Nz = supercell_dims
    
    supercell_motif = []
    for i,j,k in product(range(0,Nx),range(0,Ny),range(0,Nz)):
        for pos in motif:
            new_pos = pos + np.array([i,j,k])
            supercell_motif.append(new_pos/supercell_dims)
    
    return np.array(supercell_motif)
            

def sort_configs_by_volume(configs: list[Atoms]) -> list[Atoms]:
    
    volumes = np.array([c.get_volume() for c in configs])
    indicies = np.argsort(volumes)
    return [configs[i] for i in indicies]


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
    


            