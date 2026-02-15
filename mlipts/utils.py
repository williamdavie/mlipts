
'''
MLIPTS utilities.
'''

import numpy as np
from itertools import product
import ase

#------------------------------Reference positions------------------------------------
'''
Getting a set of reference positions that match the equilibrium configuration i.e. a high-symmetry lattice structure, 

Rephrase: We know the ground state structure of a material is fcc, but we hold a dataset with high temperature distorted supercells with non-diagonal supercell matricies. 
We want to ask, at this level of thermal expansion, what is the motif that corresponds to a perfect fcc lattice for this supercell?
    
used directly in magmom
'''


def get_scaled_reference_positions(config: ase.Atoms, equilibrium_config: ase.Atoms):
    '''
    Returns the scaled positions of a input config reformatted to align with the equilibrium configuration (high-symmetry lattice structure)
    '''

    cell = config.cell
    minimal_cell = equilibrium_config.cell


    scaling_factor = int(np.linalg.det(cell)/np.linalg.det(minimal_cell)) # intergers 1,2,3 .. etc
    

    # should be close to 1
    volume_scaling_factor = (np.linalg.det(cell)/scaling_factor)/np.linalg.det(minimal_cell)
    expanding_vectors =  (volume_scaling_factor)**(1/3) * minimal_cell

    #reference position, this is jubious but the main thing is that we select the correct element as reference
    reference_element = get_reference_element(equilibrium_config)
    reference_positions = get_reference_positions(config)
    reference_position =  reference_positions[reference_element]
    

    input_pos_ref = config.get_positions() - reference_position
    input_pos_scaled = np.zeros_like(input_pos_ref)
    
    for i,pos in enumerate(input_pos_ref):
        input_pos_scaled[i] = np.linalg.solve(expanding_vectors,pos)
        
    return input_pos_scaled



def get_reference_positions(atoms: ase.Atoms):
    
    reference_positions = {}
    dict_keys = []
    atoms_symbols = atoms.get_chemical_symbols()
    atoms_positions = atoms.get_positions()
    
    for i,element in enumerate(atoms_symbols):
        if element not in dict_keys:
            dict_keys.append(element)
            
            elem_positions = atoms_positions[[s == element for s in atoms_symbols]]
            closest_pos = elem_positions[np.argmin(np.linalg.norm(elem_positions, axis=1))]
            
            reference_positions[str(element)] = closest_pos
            
    return reference_positions


def get_reference_element(atoms: ase.Atoms):
    
    pos = atoms.positions
    idx = np.argmin(np.linalg.norm(pos, axis=1))

    closest_pos = pos[idx]
    closest_element = atoms[idx].symbol
    
    return closest_element

#-------------------------------------------------------------------------------------

def sort_configs_by_volume(configs: list[ase.Atoms]) -> list[ase.Atoms]:
    '''
    Sort a set of configurations by volume
    '''
    
    volumes = np.array([c.get_volume() for c in configs])
    indicies = np.argsort(volumes)
    return [configs[i] for i in indicies]

#-------------------------------Supercells---------------------------------------------

def match_config_to_dir(config: ase.Atoms, supercell_dict: dict) -> str:
    '''
    Given a directory, labelled by a set of lattice vectors, match the config to a directory.
    '''
    
    dirs = list(supercell_dict.keys())
    
    min_val = 100

    for i,value in enumerate(supercell_dict.values()):
        
        dif = abs(np.linalg.norm(value-np.array(config.cell)))
        print(dif)
        if dif < min_val:
            min_val = dif
            index = i
    
    if min_val > np.min(supercell_dict.values()[index])/2:
        print('Warning <!>: Did not find a matching supercell in the input dictionary')
        
    return dirs[index]
    
            

def fetch_supercell_motif(motif: np.ndarray, supercell_dims: np.ndarray):
    
    Nx,Ny,Nz = supercell_dims
    
    supercell_motif = []
    for i,j,k in product(range(0,Nx),range(0,Ny),range(0,Nz)):
        for pos in motif:
            new_pos = pos + np.array([i,j,k])
            supercell_motif.append(new_pos/supercell_dims)
    
    return np.array(supercell_motif)

#---------------------------------------------------------------------------