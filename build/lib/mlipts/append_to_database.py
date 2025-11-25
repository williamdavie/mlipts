'''

@Author William Davie

With a VASP-DFT (or other-future versions) calculation done, we now want to append positions, energies and forces to the training database.

py4vasp makes this very straightforward. 

Developer notes: 

- Formated for MACE, details on data format:  https://github.com/ilyes319/mace-tutorials

- Currently only supports a periodic lattice

- Assumes all configs have same periodic boundry conditions.


Now __main__ is vasp specific - will need to alter if other codes added in future. 

'''

import py4vasp
import numpy as np
import sys,os

def append_vasp_calculation(vasp_dir: str, database_file: str, pbc: str='T T T'):
    
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
        
        
        
    
    #'Properties=species:S:1:pos:R:3:forces:R:3 energy=-1602.093208085042 pbc="T T T"
    

if __name__ == '__main__': 
    
    vasp_calc = sys.argv[1]
    
    database_file = sys.argv[2]
    
    append_vasp_calculation(vasp_calc,database_file)
    
    print(f'Data from {vasp_calc} has been append to the database: {database_file}')