'''
@author: William Davie

File containing vasp specific functionality. Used to build many vasp calculations.

some functionality may be generalised if other codes are added 
'''

import ase
import ase.io
import numpy as np
import shutil
import py4vasp
from pathlib import Path
from mlipts import utils
from itertools import product
import ase.build


def build_vasp_calculation(vasp_base_dir: str, config: ase.Atoms, calc_name: str, outdir: str) -> str: 
    '''
    Builds a vasp calculation directory for a given atomic configuration.
    
    Parameters
    ----------
    vasp_base_dir: str
        path to base directory, should contain POTCAR, KPOINTS and INCAR.
    config: :class:`ase.Atoms` 
        atomic configuration.
    outname: stls
        name of calculation directory
    outdir: str
        output path of calculation directory.
        
    Returns
    -------
    new_calc_dir: str 
        vasp directory generated in outdir.
    '''
    
    new_calc_dir = outdir + '/' + f'{calc_name}'
    shutil.copytree(vasp_base_dir, new_calc_dir, dirs_exist_ok=True)
    ase.io.write(new_calc_dir +'/POSCAR', config, format='vasp', vasp5=True, direct=True)
    
    return new_calc_dir


def write_POSCAR_str(config: ase.Atoms) -> str:
    '''
    writes a POSCAR string given an atomic configuration.
    '''
    
    poscar = 'System\n 1.0\n'
    
    cell = np.array(config.cell)
    poscar += f' {cell[0,0]} {cell[0,1]} {cell[0,2]}\n {cell[1,0]} {cell[1,1]} {cell[1,2]}\n {cell[2,0]} {cell[2,1]} {cell[2,2]}\n'
    
    type_list = list(config.symbols)
    
    # set can be unordered so can't use set(config.symbols)
    type_labels = []
    for i in type_list:
        if i not in type_labels:
            type_labels.append(i) # only way i can see to gaurentee order?
    
    for type in type_labels:
        poscar += f' {type} '
    poscar+='\n'

    for type in type_labels:
        count = config.symbols.count(type)
        poscar += f' {count} '
    poscar+='\nCartesian\n'

    for pos in config.get_positions():
        poscar+=f'{pos[0]} {pos[1]} {pos[2]}\n'
 
    return poscar


def append_vasp_calc_to_database(database_file: str, vasp_dir: str):
    atoms = ase.io.read(f"{vasp_dir}/vasprun.xml")
    outcar_str = open(f'{vasp_dir}/OUTCAR','r').read()
    if 'aborting loop EDIFF was not reached (unconverged)' in outcar_str:
        print('Self consistency failed, not saving data.')
        return None
    ase.io.write(database_file, atoms, format="extxyz", append=True)
    return None



def fetch_configs_vasp(calc_dirs: list[str]) -> list[ase.Atoms]:
    '''
    from a set of directories containing vasp in files, read configs.
    '''
    configs = []
    for dir in calc_dirs:
        if havePOSCAR(dir):
            atoms = ase.io.read(f'{dir}/POSCAR')
            configs.append(atoms)
        else:
            print(f'No POSCAR found in directory: {dir}')
            pass
        
    return configs 
            
        


'''
Want some native way of editing vasp calculations. Namely (for my work) increasing magmom for supercells. 
'''

#-------------------set ICHARG across database------------------

def set_icharg(value: int, vasp_calc_dir: str):
    
    incar_lines = open(f'{vasp_calc_dir}/INCAR','r').readlines()
    found = False
    for i,line in enumerate(incar_lines):
        if 'ICHARG' in line:
            incar_lines[i] = f'ICHARG = {value}\n'
            found = True
    if not found:
        incar_lines.append('\n')
        incar_lines.append(f'ICHARG = {value}\n')
    
    with open(f'{vasp_calc_dir}/INCAR','w') as f:
        new_file_str = "".join(incar_lines)
        f.write(new_file_str)
        
    return None



#-------------------MAGMOM for large databases------------------


def set_magmom(magmom_motif_config: ase.Atoms, vasp_calc_dirs: str='./QM_calculations') -> None:
    '''
    Given a set of vasp calculation directories, the supercell size, a motif and the magnet moments for the motif, POSCAR is used to set the MAGMOM string. 
    Allowing the user to access magnetically ordered states for larger supercells.
    
    This is a solid specific functionality where 
    
    Parameters
    ----------
    supercell_size :class:`np.ndarray` 
        3D array defining supercell size
    motif: :class:`np.ndarray` 
        motif of a relaxed solid structure. 
    magmom_motif: :class:`np.ndarray` 
        magnetic moments of the motif, order of magmom_motif must equal the order of motif. i.e. the magnetic moment of atom located at motif[i] is magmom_motif[i].
        
    Returns
    -------
    None : None
        edits INCAR files in call sub directories. 
    '''
    
    path = Path(vasp_calc_dirs)
    subdirs = [p for p in path.iterdir() if p.is_dir()]
    for vasp_calc in subdirs:
        if haveINCAR(str(vasp_calc)) and havePOSCAR(str(vasp_calc)):
            set_magmom_one_directory(magmom_motif_config,vasp_calc)
        else:
            pass
    
    print(f'Magnetic Moments updated in all vasp sub directories of {vasp_calc_dirs}')
    
    return None
    
    
def set_magmom_one_directory(magmom_motif_config: ase.Atoms, vasp_calc_dir: str) -> None:
    '''
    Called on each directory by set_magmom
    '''
    # This function is quite brute force and is oppitunity to optimize.
    
    # read data
    minimal_cell = np.array(magmom_motif_config.cell)
    try:
        motif_magmoms = magmom_motif_config.get_initial_magnetic_moments()
        motif_scaled_pos = magmom_motif_config.get_scaled_positions()
        motif_elements = magmom_motif_config.get_chemical_symbols()
    except:
        raise ValueError('magnetic moment motif input file must contain positions and initial moments.')

    input_atoms = ase.io.read(f'{vasp_calc_dir}/POSCAR')
    cell = np.array(input_atoms.cell)
    
    # create a high symmetry configuration for an expanded cell
    expanded_magmom_motif_config = magmom_motif_config.copy()

    scaling_factor = round(np.linalg.det(cell)/np.linalg.det(minimal_cell))
    volume_scaling_factor = (np.linalg.det(cell)/scaling_factor)/np.linalg.det(minimal_cell)
    expanded_cell =  (volume_scaling_factor)**(1/3) * minimal_cell
    
    expanded_magmom_motif_config.set_cell(expanded_cell)
    expanded_magmom_motif_config.set_scaled_positions(motif_scaled_pos)
    expanded_magmom_motif_config.set_initial_magnetic_moments(motif_magmoms)
    
    # now expand this cell to the input supercell.
    
    S = utils.get_supercell_matrix(cell,expanded_cell)
    
    if np.allclose(S, np.round(S), atol=1e-4):
        supercell_magmom_motif = ase.build.make_supercell(expanded_magmom_motif_config,np.round(S))
    else:
        raise ValueError('cannot set magmom for configurations with non-interger supercell matricies, maybe need to use a conversion before using as input data.')
    
    final_cell = supercell_magmom_motif.cell
    final_elements =supercell_magmom_motif.get_chemical_symbols()
    final_magmoms =supercell_magmom_motif.get_initial_magnetic_moments()
    
    possible_vectors = []
    # by expanding range to (-1,1), variations of wrapped co-ords outputed by the MD calculation are dealt with.
    for i,j,k in product(range(-1,2),range(-1,2),range(-1,2)):
        possible_vectors.append(np.array([i,j,k]))
    expected_positions = [] # expected for a relaxed lattice
    expected_mag_moments = []
    expected_elements = []
    for vecs in possible_vectors:
        for l,equilibrium_pos in enumerate(supercell_magmom_motif.get_positions()):
            pos = equilibrium_pos + vecs @ final_cell
            expected_positions.append((pos))
            expected_elements.append(final_elements[l])
            expected_mag_moments.append(final_magmoms[l]) # set corresponding magmom
   
    # Now match magmom according to minimum distance and element.
    
    A = input_atoms.positions
    B = np.array(expected_positions)
    all_diff = B[None,:,:] - A[:, None, :]
    all_dist2 = np.sum(all_diff**2, axis=2)  
    
    magmoms = np.zeros((len(input_atoms),3)) 
    
    for i,atom in enumerate(input_atoms):
        symbol = atom.symbol
        element_indices = [i for i, val in enumerate(expected_elements) if val == symbol]
        B_this_element = B[np.array(element_indices), :]
        diff = B_this_element - A[i]
        dist2 = np.sum(diff**2, axis=1)  
        closest_index = np.argmin(dist2)
        true_index = np.where(dist2[closest_index]==all_dist2)[1][0]
        magmoms[i] = expected_mag_moments[true_index]
    
    # with magmom known can define magmom str
    magmom_str = 'MAGMOM = '
    for i, pos in enumerate(input_atoms.positions):
        mx,my,mz = magmoms[i][0:3]
        magmom_str += f'{mx} {my} {mz} '
        
    input_atoms.set_initial_magnetic_moments(magmoms)
    writeMAGMOM(f'{vasp_calc_dir}/INCAR',new_magmom_str=magmom_str)
    
    return None
            
    
def writeMAGMOM(incar: str, new_magmom_str: str) -> None:
    '''
    given the path to an INCAR file, writes or updates the MAGMOM string
    '''
    
    incar_lines = open(incar,'r').readlines()
    found = False
    for i,line in enumerate(incar_lines):
        if 'MAGMOM' in line:
            incar_lines[i] = new_magmom_str + '\n'
            found = True
    if not found:
        incar_lines.append('\n')
        incar_lines.append(new_magmom_str + '\n')
        
    with open(incar,'w') as f:
        new_file_str = "".join(incar_lines)
        f.write(new_file_str)
        
    return None



#-------------------KPOINTS for large databases------------------

def set_kpoints(kspacing: float, grid_type: str='Gamma',vasp_calc_dirs: str='./QM_calculations'):
    
    path = Path(vasp_calc_dirs)
    subdirs = [p for p in path.iterdir() if p.is_dir()]
    for vasp_calc in subdirs:
        if havePOSCAR(str(vasp_calc)):
            set_kpoints_one_directory(str(vasp_calc),kspacing,grid_type)
        else:
            pass

    return None
    
    

def set_kpoints_one_directory(vasp_calc_dir: str, kspacing: float, grid_type: str='Gamma'):
    '''

    :param kspacing: units: 2pi/A

    '''
    
    input_atoms = ase.io.read(f'{vasp_calc_dir}/POSCAR')
    cell = np.array(input_atoms.cell)
    rcp_lattice_vectors = reciprocal_lattice_vectors(cell)
    
    kpoints = np.zeros(3)
    for i in range(3):
        # There are some subtlties 
        kpoints[i] = max(1,round(np.linalg.norm(rcp_lattice_vectors[i])/(2*np.pi*kspacing)))
    
    
    with open(f'{vasp_calc_dir}/KPOINTS','w') as f:
        
        f.write(f'K-Spacing Value to Generate K-Mesh: {kspacing:.3f}\n')
        f.write('0\n')
        f.write(f'{grid_type}\n')
        f.write(f'{kpoints[0]:.0f} {kpoints[1]:.0f} {kpoints[2]:.0f}\n')
        f.write('0.0 0.0 0.0\n')
        
    return None



def reciprocal_lattice_vectors(lattice_vectors: np.ndarray):
    '''
    Returns the reciprocal lattice vectors of 3x3 lattice vectors.
    '''
    
    rcp_lattice_vectors = np.zeros_like(lattice_vectors)
    V = np.dot(lattice_vectors[0], np.cross(lattice_vectors[1],lattice_vectors[2]))
    
    pertubations = [(0,1,2),(1,2,0),(2,0,1)]
    for p in pertubations:
        rcp_lattice_vectors[p[0]] = 2*np.pi/V * (np.cross(lattice_vectors[p[1]],lattice_vectors[p[2]]))
        
    return rcp_lattice_vectors

#-------------------ANY INCAR PARAM for large databases------------------

# want a way to set an INCAR parameter given a cell type condition, e.g. number of atoms.

    
    
def haveINCAR(dir: str):
    '''
    checks if a directory contains INCAR
    '''
    path = Path(dir)
    files = [str(p.name) for p in path.iterdir()]
    if 'INCAR' in files:
        return True
    else:
        return False
    
def havePOSCAR(dir: str):
    '''
    checks if a directory contains POSCAR
    '''
    path = Path(dir)
    files = [str(p.name) for p in path.iterdir()]
    if 'POSCAR' in files:
        return True
    else:
        return False
    
def haveKPOINTS(dir: str):
    '''
    checks if a directory contains KPOINTS
    '''
    path = Path(dir)
    files = [str(p.name) for p in path.iterdir()]
    if 'KPOINTS' in files:
        return True
    else:
        return False
        
        
    

    
    

    