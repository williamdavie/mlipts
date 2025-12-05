'''
@author: William Davie

File containing vasp specific functionality. Used to build many vasp calculations.
'''

from ase import Atoms
from ase.io import read, write
import numpy as np
import shutil


def build_vasp_calculation(vasp_base_dir: str, config: Atoms, calc_name: str, outdir: str) -> None: 
    '''
    Builds a vasp calculation directory for a given atomic configuration.
    
    Parameters
    ----------
    vasp_base_dir: str
        path to base directory, should contain POTCAR, KPOINTS and INCAR.
    config: :class:`ase.Atoms` 
        atomic configuration.
    outname: str
        name of calculation directory
    outdir: str
        output path of calculation directory.
        
    Returns
    -------
    None: None
        vasp directory generated in outdir.
    '''
    
    poscar = write_POSCAR_str(config)
    new_calc_dir = outdir + '/' + f'{calc_name}'
    shutil.copytree(vasp_base_dir, new_calc_dir, dirs_exist_ok=True)
            
    with open(new_calc_dir +'/POSCAR','w') as f:
                f.write(poscar)
    
    return None


def write_POSCAR_str(config: Atoms) -> str:
    '''
    writes a POSCAR string given an atomic configuration.
    '''
    
    poscar = 'System\n 1.0\n'
    
    cell = np.array(config.cell)
    poscar += f' {cell[0,0]} {cell[0,1]} {cell[0,2]}\n {cell[1,0]} {cell[1,1]} {cell[1,2]}\n {cell[2,0]} {cell[2,1]} {cell[2,2]}\n'
    
    type_labels = config.symbols.species()
    for type in type_labels:
        poscar += f' {type} '
    poscar+='\n'

    for type in type_labels:
        count = config.symbols.count(type)
        poscar += f' {count} '
    poscar+='\nDirect\n'
    

    for pos in config.positions:
        poscar+=f'{pos[0]} {pos[1]} {pos[2]}\n'
 
    return poscar
             