'''
@author: William Davie

File containing vasp specific functionality. Used to build many vasp calculations.

some functionality may be generalised if other codes are added 
'''

from ase import Atoms
from ase.io import read, write
import numpy as np
import shutil
import py4vasp


def build_vasp_calculation(vasp_base_dir: str, config: Atoms, calc_name: str, outdir: str) -> str: 
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
    new_calc_dir: str 
        vasp directory generated in outdir.
    '''
    
    poscar = write_POSCAR_str(config)
    new_calc_dir = outdir + '/' + f'{calc_name}'
    shutil.copytree(vasp_base_dir, new_calc_dir, dirs_exist_ok=True)
            
    with open(new_calc_dir +'/POSCAR','w') as f:
                f.write(poscar)
    
    return new_calc_dir


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
             
        
def append_vasp_calc_to_database(database_file: str, vasp_dir: str, pbc: str='T T T'):
    
    final_config = ''

    try:
        calc = py4vasp.Calculation.from_path(vasp_dir)
        # fetch energy
        energy_toten = calc.energy.to_dict()
        eval = energy_toten['free energy    TOTEN']

        # fetch force and structure data
        force_data =  calc.force.to_dict()
        elements = force_data['structure']['elements']
        num_atoms = len(elements)

        l_vecs = force_data['structure']['lattice_vectors']
        lattice_str = f'Lattice="{l_vecs[0,0]} {l_vecs[0,1]} {l_vecs[0,2]} {l_vecs[1,0]} {l_vecs[1,1]} {l_vecs[1,2]} {l_vecs[2,0]} {l_vecs[2,1]} {l_vecs[2,2]}"'
        atomic_positions = force_data['structure']['positions']
    
        forces = force_data['forces']
        final_config += f'{num_atoms}\n'
        final_config += f'{lattice_str} Properties=species:S:1:pos:R:3:forces_xtb:R:3 energy_xtb={eval} pbc="{pbc}"\n'
    
        for i in range(0,num_atoms):
            # convert to cartiesian
            species = elements[i]
            position = atomic_positions[i][0] * l_vecs[0] + atomic_positions[i][1] * l_vecs[1] + atomic_positions[i][2] * l_vecs[2]
            force = forces[i] # possibly a conversion required.
            final_config += f'{species} {position[0]} {position[1]} {position[2]} {force[0]} {force[1]} {force[2]}\n'
        
        with open(database_file,'a') as f:
            f.write(final_config)
            
    except Exception as e:
        print(f'The calculation under {vasp_dir} did not failed or did not finish.')
        print(e)
        return None