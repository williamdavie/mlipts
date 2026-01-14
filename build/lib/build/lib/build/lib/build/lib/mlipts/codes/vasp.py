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

from itertools import product


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
    
    poscar = write_POSCAR_str(config)
    new_calc_dir = outdir + '/' + f'{calc_name}'
    shutil.copytree(vasp_base_dir, new_calc_dir, dirs_exist_ok=True)
            
    with open(new_calc_dir +'/POSCAR','w') as f:
                f.write(poscar)
    
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
    poscar+='\nDirect\n'

    for pos in config.get_scaled_positions():
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


def set_magmom(supercell_size: np.ndarray, 
               motif: np.ndarray, magmom_motif: np.ndarray, 
               vasp_calc_dirs: str='./QM_calculations') -> None:
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
            set_magmom_one_directory(supercell_size,motif,magmom_motif,vasp_calc)
        else:
            pass
    
    print(f'Magnetic Moments updated in all vasp sub directories of {vasp_calc_dirs}')
    
    return None
    
    
def set_magmom_one_directory(supercell_size: np.ndarray,
                             motif: np.ndarray, magmom_motif: np.ndarray,
                             vasp_calc_dir: str) -> None:
    '''
    Called on each directory by set_magmom
    '''
    # This function is quite brute force and is oppitunity to optimize.
    
    # define all possible positions
    atoms = ase.io.read(f'{vasp_calc_dir}/POSCAR')
    basis_vectors = np.array(atoms.cell)/supercell_size

    Nx,Ny,Nz = supercell_size[0:3]
    possible_vectors = []
    # by expanding range to (-1,N+1, variations of wrapped co-ords outputed by the MD calculation. 
    for i,j,k in product(range(-1,Nx+1),range(-1,Ny+1),range(-1,Nz+1)):
        possible_vectors.append(np.array([i,j,k]))
    expected_positions = [] # expected for a relaxed lattice
    mag_moments = []
    for vecs in possible_vectors:
        for i,motif_pos in enumerate(motif):
            pos = (motif_pos + vecs)
            pos_cart = pos[0] * basis_vectors[0] + pos[1] * basis_vectors[1] + pos[2] * basis_vectors[2]
            expected_positions.append((pos_cart))
            mag_moments.append(magmom_motif[i]) # set corresponding magmom
            
    # find the positions in POSCAR corresponding to positions in motif
    A = atoms.positions
    B = np.array(expected_positions)
    diff = A[:, None, :] - B[None, :, :]  
    dist2 = np.sum(diff**2, axis=2)       
    closest_indices = np.argmin(dist2, axis=1)
    magmom_reordered = np.array(mag_moments)[closest_indices]
    
    # define magmom str
    magmom_str = 'MAGMOM = '
    for i, pos in enumerate(atoms.positions):
        mx,my,mz = magmom_reordered[i][0:3]
        magmom_str += f'{mx} {my} {mz} '
        
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
        
    

    
    

    