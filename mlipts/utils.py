
'''
MLIPTS utilities.
'''

import numpy as np
from itertools import product
import ase
import scipy.spatial

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


    
def return_motif_config(config: ase.Atoms, equilibrium_config: ase.Atoms) -> ase.Atoms:
    '''
    Some positions may be wrapped to larger cell sizes. 
    
    <!> can generalise with set_magmom() in future.
    '''
    motif_config = config.copy()
    
    cell = config.cell
    minimal_cell = equilibrium_config.cell
    # create a high symmetry configuration for an expanded cell
    expanded_config = equilibrium_config.copy()

    scaling_factor = round(np.linalg.det(cell)/np.linalg.det(minimal_cell))
    volume_scaling_factor = (np.linalg.det(cell)/scaling_factor)/np.linalg.det(minimal_cell)
    expanded_cell =  (volume_scaling_factor)**(1/3) * minimal_cell
    
    expanded_config.set_cell(expanded_cell)
    expanded_config.set_scaled_positions(equilibrium_config.get_scaled_positions())
    
    # now expand this cell to the input supercell.
    
    S = get_supercell_matrix(cell,expanded_cell)
    
    if np.allclose(S, np.round(S), atol=1e-4):
        supercell_motif = ase.build.make_supercell(expanded_config,np.round(S))
    else:
        raise ValueError('cannot set magmom for configurations with non-interger supercell matricies, maybe need to use a conversion before using as input data.')

    positions_extended = []
    for i,j,k in product(range(-1,2),range(-1,2),range(-1,2)):
        for l,pos in enumerate(supercell_motif.get_positions()):
            new_pos = pos + np.array([i,j,k]) @ supercell_motif.cell
            #pos_cart = pos[0] * lattice_vectors[0] + pos[1] * lattice_vectors[1] + pos[2] * lattice_vectors[2]
            positions_extended.append(new_pos)
            
    # notice this is very similar to that used in mlipts.codes.vasp.set_magmom, could be generalized. 
    A = config.get_positions()
    B = np.array(positions_extended)
    diff =  B[None, :, :]  - A[:,None,:]
    dist2 = np.sum(diff**2, axis=2)       
    closest_indices = np.argmin(dist2, axis=1)
    positions_new = B[closest_indices]
    
    motif_config.set_positions(positions_new)

    return motif_config


#-------------------------------Defects---------------------------------------------


def generate_defect(config: ase.Atoms, targets: dict[str, int], defect_type: str='schottkey') -> ase.Atoms:
    '''
    Fenerates a config that contains a defect.
    
    Parameters
    ----------
    
    config: ase.Atoms
        an initial configuration of atoms
    targets: dict[str, int]
        dictionary defining the elements involved with defect. Format: 'element': count
    defect_type: str
        defines the defect type 'schottkey' or 'frenkel'. 
    
    Returns
    -------

    config: ase.Atoms
        defected atomic configuration.
    '''

    if defect_type == 'schottkey':
        config = generate_schottkey_defect(config,targets)
    elif defect_type == 'frenkel':
        config = generate_frenkel_defect(config,targets)
    
    
    return config


def generate_schottkey_defect(config: ase.Atoms, targets: dict[str, int]):
    '''generate schottkey defect'''
    
    symbols = np.array(config.get_chemical_symbols(),dtype='U10')
    
    # find smallest atom count (should have highest charge).
    minimum_count = min(targets.values()) 
    minimum_key = [key for key in targets if targets[key]==minimum_count][0] # if more than one value just select the first.
    
    # choose largest elements
    old_config = config.copy()
    symbols = np.array(config.get_chemical_symbols(),dtype='U10')
    removed_largest = []
    min_key_indices = np.where(symbols == minimum_key)[0]
    for l in range(minimum_count):
        
        del config[min_key_indices[l]]
        removed_largest.append(min_key_indices[l])

    # remove nearest neighbours of highest charge element.
    surrounding_elements = [key for key in targets if targets[key]!=minimum_count]
    for j in surrounding_elements:
        to_remove = find_nearest_neighbours(old_config,removed_largest,Nneighbours=targets[j],target_element=j)
        for k in to_remove:
            del config[k]
            
    return config

    
def generate_frenkel_defect(config: ase.Atoms, targets: dict[str, int]):
    '''generate frenkel defect'''
    
    symbols = np.array(config.get_chemical_symbols(),dtype='U10')
    
    for target_element in targets.keys():
        atom_indices = np.where(symbols == target_element)[0]
        count = targets[target_element]
        interstitial_sites = find_interstitial_sites(config, count)
        for i in range(count):
            config.positions[atom_indices[i]] = interstitial_sites[i]
        
    return config


def find_nearest_neighbours(config: ase.Atoms, central_atoms: list[int], Nneighbours: int, target_element: str='O'):
    '''Finds nearest N neighbours of central_atoms (indicies of the config.)'''
    
    symbols = np.array(config.get_chemical_symbols(),dtype='U10')
    atom_indices = np.where(symbols == target_element)[0]
    
    atom_positions = config.get_positions()[atom_indices]
    
    original_indices = []
    distances = []
    for i in central_atoms:
        target_position = config.positions[i]
        distances.extend(np.linalg.norm(atom_positions - target_position,axis=1))
        original_indices.extend([i for i in range(len(atom_positions))])

    min_distances = np.array(original_indices)[np.argsort(distances)]

    atom_indices = atom_indices[min_distances]
    
    return atom_indices[0:Nneighbours] 
    
    
def find_interstitial_sites(atoms: ase.Atoms, count: int, grid_density: int=10):
    '''Basic function to find an interstitial site.'''
    
    pos = atoms.get_positions()
    cell = atoms.get_cell()
    
    x = np.linspace(0, 0.5, grid_density)
    grid_points_frac = np.array(np.meshgrid(x, x, x)).T.reshape(-1, 3)
    grid_points_cart = np.dot(grid_points_frac, cell)
    
    tree = scipy.spatial.KDTree(pos)
    
    distances, _ = tree.query(grid_points_cart)
    
    best_idx = np.argmax(distances)
    best_indices = np.argsort(distances)[::-1]

    best_sites = grid_points_cart[best_indices[0:count]]

    return best_sites


    