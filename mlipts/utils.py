import numpy as np
from itertools import product
import ase

def fetch_supercell_motif(motif: np.ndarray, supercell_dims: np.ndarray):
    
    Nx,Ny,Nz = supercell_dims
    
    supercell_motif = []
    for i,j,k in product(range(0,Nx),range(0,Ny),range(0,Nz)):
        for pos in motif:
            new_pos = pos + np.array([i,j,k])
            supercell_motif.append(new_pos/supercell_dims)
    
    return np.array(supercell_motif)
            

def sort_configs_by_volume(configs: list[ase.Atoms]) -> list[ase.Atoms]:
    
    volumes = np.array([c.get_volume() for c in configs])
    indicies = np.argsort(volumes)
    return [configs[i] for i in indicies]


def match_config_to_dir(config: ase.Atoms, supercell_dict: dict) -> str:
    '''
    Given a directory, labelled by a supercell matrix and an atomic config, match the config to a directory.
    '''
    
    keys = list(supercell_dict.keys())
    if 'a' not in keys:
        raise ValueError('Tried to a config to the corresponding directory but lattice parameter (a) was not provided')

    index = 0
    min_val = supercell_dict['a'] * 2

    for i,value in enumerate(supercell_dict.values()):
        
        if keys[i] == 'a':
            continue
        
        print(value*supercell_dict['a'])
        print(np.array(config.cell))
        dif = abs(np.linalg.norm(value*supercell_dict['a']-np.array(config.cell)))
        print(dif)
        if dif < min_val:
            min_val = dif
            index = i
    
    if min_val > supercell_dict['a']/2:
        print('Warning <!>: Did not find a matching supercell in the input dictionary')
        
    return list(supercell_dict.keys())[index]


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
    


            