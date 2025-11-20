'''
@Author William Davie

Here we use a dump output from LAMMPS and construct VASP inputs (POSCAR files)

- The user will input a lammps directory

- Within that directory a number of sub directories will be created for vasp calculations.

- Currently requires user to also provide a 'vasp base directory', this includes POTCAR KPOINTS and INCAR

Development notes:

- Hopefully will be generalized for any MD code, but currently only support LAMMPS

- writing POSCAR currently relies on an orthorhombic simulation box

'''

from pathlib import Path
import numpy as np
from typing import List
import shutil
import os

class constructVASP():
    
    def __init__(self, vasp_base_dir: str, md_dir: str):
        
        self.vasp_base_dir = vasp_base_dir
        
        # check VASP base dir:
        
        vasp_files = os.listdir(self.vasp_base_dir)
        
        required = {'POTCAR', 'KPOINTS', 'INCAR'}
        assert required.issubset(vasp_files), 'Vasp base directory must contain POTCAR, KPOINTS and INCAR'
        
        # molecular dynamics directory
        self.md_dir = md_dir
        
        # check MD dir
    
        # this is LAMMPS specific but easy to add other options. 
        self.md_output = list(Path(self.md_dir).glob('md.*'))
        
        assert self.md_output, '''Error: MD simulation directory must contain an output file formated as "md.*"
        - Ensure you have run the relevant MD calculations. 
        '''
        
        assert len(self.md_output) == 1, 'Error: MD simulation directory has multiple files named "md.* '
    
        self.md_output = open(self.md_output[0],'r').read()
        
        
    
    def read_pos_LAMMPS(self):
        '''
        From file md.* reads output as specified in LAMMPS. 
        '''
        
        split_lines = self.md_output.splitlines()
        
        self.timesteps = [] # snapshot time steps.
        self.all_atomic_positions = [] 
        self.all_lattice_vectors = []
        
        num_snapshots = self.md_output.count('ITEM: TIMESTEP')
        self.num_atoms = int(split_lines[split_lines.index('ITEM: NUMBER OF ATOMS') + 1])
        
        #atom_information_string = split_lines[split_lines.index('')]
        
        lattice_vectors = np.zeros((3,3))
         
        for i,line in enumerate(split_lines):
            
            if 'ITEM: TIMESTEP' in line:
                
                self.timesteps.append(split_lines[i+1])

            if 'ITEM: BOX BOUNDS' in line:
                     
                aline = split_lines[i+1].split()
                bline = split_lines[i+2].split()
                cline = split_lines[i+3].split()
                
                # hi - lo
                lattice_vectors[0,0] = float(aline[1]) - float(aline[0])
                lattice_vectors[1,1] = float(bline[1]) - float(bline[0])
                lattice_vectors[2,2] = float(cline[1]) - float(cline[0])
                
    
                self.all_lattice_vectors.append(lattice_vectors.copy())
        
            
            if 'ITEM: ATOMS' in line:
                
                atomic_positions = np.zeros(shape=(self.num_atoms,3))
                
                types = np.zeros(shape=(self.num_atoms))
                
                line_split = line.split()
                
                args = line_split[2:]
                
                id_index = args.index('id')
                type_index = args.index('type')
                x_index = args.index('x')
                y_index = args.index('y')
                z_index = args.index('z')
                
                for atom in range(0,self.num_atoms):
            
                    atom_line_split = split_lines[i+atom+1].split()

                    types[atom] = int(atom_line_split[type_index]) 
                    
                    atomic_positions[atom,0] = float(atom_line_split[x_index]) / lattice_vectors[0,0]
                    atomic_positions[atom,1] = float(atom_line_split[y_index]) / lattice_vectors[1,1]
                    atomic_positions[atom,2] = float(atom_line_split[z_index]) / lattice_vectors[2,2]
                    
                # as required by vasp, must be type sorted. 
                types_sorted = np.argsort(types)

                self.types,self.type_counts = np.unique(types,return_counts=True)
        
                self.num_types = np.max(self.types)
                
                self.all_atomic_positions.append(atomic_positions[types_sorted])
                
    
    def write_POSCAR(self, lattice_vectors: np.ndarray, atomic_positions: np.ndarray, type_labels: List[str]) -> str:
        '''
        writes a POSCAR string
        '''
        
        assert len(type_labels)==self.num_types, 'Error, found more types in md.* trajectory file than number of type labels provided'
        
        poscar = 'System\n 1.0\n'
        # again requires orthohomic box
        poscar += f' {lattice_vectors[0,0]} 0.0 0.0\n 0.0 {lattice_vectors[1,1]} 0.0\n 0.0 0.0 {lattice_vectors[2,2]}\n'
    
        for type in type_labels:
            poscar += f' {type} '
        poscar+='\n'

        for i in self.type_counts:
            poscar += f' {i} '
        poscar+='\nDirect\n'
        
        self.num_atoms
        for i in range(0,int(self.num_atoms)):
            poscar+=f'{atomic_positions[i,0]} {atomic_positions[i,1]} {atomic_positions[i,2]}\n'
            
        
        return poscar
    

    def build_vasp_dirs(self, type_labels: str):
        
        self.new_vasp_calc_dirs = []
        
        for i,step in enumerate(self.timesteps):
            
            poscar = self.write_POSCAR(self.all_lattice_vectors[i],self.all_atomic_positions[i],type_labels)
            
            new_vasp_calc_dir = self.md_dir + '/' + f'vasp_step_{step}'
            
            self.new_vasp_calc_dirs.append(new_vasp_calc_dir)
            
            shutil.copytree(self.vasp_base_dir, new_vasp_calc_dir, dirs_exist_ok=True)
            
            with open(new_vasp_calc_dir +'/POSCAR','w') as f:
                f.write(poscar)
             

def vasp_for_lammps(vasp_base_dir: str,lammps_directories: list[str], atom_type_labels: list[str]):
    
    for lammps_calc in lammps_directories:
        
        build = constructVASP(lammps_calc)
        build.read_pos_LAMMPS()
        build.build_vasp_dirs(vasp_base_dir,atom_type_labels)
        
        
    
    

if __name__ == '__main__':
    test = constructVASP('./base_vasp','./Examples/example_T_100')
    test.read_pos_LAMMPS()
    test.build_vasp_dirs(['O','Pu'])
    
        
        