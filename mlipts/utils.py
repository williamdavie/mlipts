
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

#------------------------------Retrieving supercell matricies--------------------------------------

def get_supercell_matrix(cell,minimal_cell):
    '''
    Get the transformation between the original high symmetry primitive lattice vectors and the input cell.
    '''

    scaling_factor = round(np.linalg.det(cell)/np.linalg.det(minimal_cell)) # intergers 1,2,3 .. etc
    
    # should be close to 1
    volume_scaling_factor = (np.linalg.det(cell)/scaling_factor)/np.linalg.det(minimal_cell)
    
    expanded_minimal_cell =  (volume_scaling_factor)**(1/3) * minimal_cell
    
    # H' = SH, solve for S
    S = np.linalg.solve(expanded_minimal_cell.T,cell.T).T
    
    return S


def retrieve_standard_supercell(cell,minimal_cell):
    '''
    Given a rotated non-diagonal supercell and the primitive lattice vectors. Calculate the interger supercell matrix.
    
    This is achieved by mapping the angles and norm ratios between vectors for all possible supercell matricies with the expected determinant. 
    '''
    
    supercell_det = round(np.linalg.det(cell)/np.linalg.det(minimal_cell))
    
    cell_parameters = get_magnitudes_and_angles(cell)
    
    possible_supercells = get_possible_supercells(supercell_det)
    
    for supercell in possible_supercells:
        
        supercell_parameters = get_magnitudes_and_angles(supercell)

        
        if np.allclose(supercell_parameters, cell_parameters, atol=1e-4):
            
            parameter = np.linalg.norm(cell)/np.linalg.norm(supercell)
           
            result = parameter * supercell
        
            return result
        

def get_possible_supercells(det: int):
    '''
    given the desired matrix determinant, return all possible upper triangular HNF supercell matricies
    '''
       
    # fetch all diagonal elements             
    diagonal_elements = []
    for i in range(1, det + 1):
        if det % i == 0:
            remaining = det // i
            for j in range(1, remaining + 1):
                if remaining % j == 0:
                    k = remaining // j
                    diagonal_elements.append((i, j, k))
              

    # define all given a diagonal element
    supercells = []
    for diag in diagonal_elements:
        S00,S11,S22 = diag
        for n in range(0,diag[1]):
            S01 = n
            for p,q in product(range(0,diag[2]),range(0,diag[2])):
                S02,S12 = p,q
                
                S = np.array([[S00,S01,S02],
                              [0,S11,S12],
                              [0,0,S22]])
                
                supercells.append(S)
                
    return supercells
                
                
def get_magnitudes_and_angles(cell):

    alpha = vectors_angle(cell[1],cell[2])
    beta = vectors_angle(cell[0],cell[2])
    gamma = vectors_angle(cell[0],cell[1])
    
    mag_array = np.array([np.linalg.norm(cell[0]),
                            np.linalg.norm(cell[1]),
                            np.linalg.norm(cell[2])])
    mag_array /= np.min(mag_array)

    return np.array([alpha,beta,gamma,mag_array[0],mag_array[1],mag_array[2]])


def vectors_angle(vec1,vec2):
    
    cosAngle = np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))
    
    return np.arccos(cosAngle)


        