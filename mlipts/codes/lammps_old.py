'''

@Author: William Davie

Automatically set's up a LAMMPS calculation, snapshots of this calculation are used for DFT caclulations. 

The user provides a 'base directory' including:

    - The input file
        - this should include the potential specification, often this will refer to a file elsewhere and a minimization task. 
        - (Option) A mlipts specific variable marked with £
        
    - The .dat file containing atomic positions
    
    - Note that all contents of the base directory will be copied. 
    
The user also provides a pytrain 'variational input file', this file includes additional lines of a LAMMPS input file to loop through. 

Notes:

- lmps shorthand for LAMMPS
- if multiple variables are specified it will generate N_var1 * N_var2 directories.

'''

import os, sys
import shutil
from pathlib import Path
import numpy as np
import re
from itertools import product
from ase import Atoms


def build_lammps_calculations(base_dir: str) -> list[str]:
    '''
    Generates a set of directories containing input files for lammps calculations
    '''
    
    build = buildLAMMPS(f'{base_dir}')
    build.fetch_variables() 
    build.set_variables()
    build.generate_outputs()
    return build.new_dirs


def write_lammps_bash_command(lammps_dirs: list[str], lammps_cmd_line: str, output_directory: str ='.'):
    '''
    Generates the command line to run a set of lammps directories.
    '''

    dirs = ''
    for dir in lammps_dirs:
        dirs += f'{dir} '
        
        input_path = list(Path(dir).glob('in.*'))[0]
        input_name = input_path.name
        
            
        output_name = 'out.' + input_name.split(".")[1]
        
        result = f'directories="{dirs}"\n'
        
        result += f'cd {output_directory}\n'
        
        result+=f'''for i in $directories; 
do 
cd $i
    {lammps_cmd_line} -i {input_name} -l {output_name}
cd ..
done\n'''
        result += 'cd ..'

        
        return result
    

def read_lammps_output(output_dir: str, atom_types: list[str], pbc: bool=True) -> list[Atoms]:
    '''
    Read lammps output.
    '''
    
    output_dir = list(Path(output_dir).glob('md.*'))
        
    
    md_output = open(output_dir[0],'r').read()
    
    
    split_lines = md_output.splitlines()
    
    timesteps = [] # snapshot time steps.
    all_atomic_positions = [] 
    all_lattice_vectors = []
    type_labels = []
    
    num_snapshots = md_output.count('ITEM: TIMESTEP')
    num_atoms = int(split_lines[split_lines.index('ITEM: NUMBER OF ATOMS') + 1])
    
    #atom_information_string = split_lines[split_lines.index('')]
    
    lattice_vectors = np.zeros((3,3))
        
    for i,line in enumerate(split_lines):
        
        if 'ITEM: TIMESTEP' in line:
            
            timesteps.append(split_lines[i+1])

        if 'ITEM: BOX BOUNDS' in line:
                    
            aline = split_lines[i+1].split()
            bline = split_lines[i+2].split()
            cline = split_lines[i+3].split()
            
            # hi - lo
            lattice_vectors[0,0] = float(aline[1]) - float(aline[0])
            lattice_vectors[1,1] = float(bline[1]) - float(bline[0])
            lattice_vectors[2,2] = float(cline[1]) - float(cline[0])
            

            all_lattice_vectors.append(lattice_vectors.copy())
    
        
        if 'ITEM: ATOMS' in line:
            
            atomic_positions = np.zeros(shape=(num_atoms,3))
            
            types = np.zeros(shape=(num_atoms))
            
            line_split = line.split()
            
            args = line_split[2:]
            
            id_index = args.index('id')
            type_index = args.index('type')
            x_index = args.index('x')
            y_index = args.index('y')
            z_index = args.index('z')
            
            for atom in range(0,num_atoms):
        
                atom_line_split = split_lines[i+atom+1].split()

                types[atom] = int(atom_line_split[type_index]) 
                
                atomic_positions[atom,0] = float(atom_line_split[x_index]) / lattice_vectors[0,0]
                atomic_positions[atom,1] = float(atom_line_split[y_index]) / lattice_vectors[1,1]
                atomic_positions[atom,2] = float(atom_line_split[z_index]) / lattice_vectors[2,2]
                
            # as required by vasp, must be type sorted. 
            types_sorted = np.argsort(types)

            types,type_counts = np.unique(types,return_counts=True)
        
            num_types = np.max(types)
            
            all_atomic_positions.append(atomic_positions[types_sorted])

    
    for i in types:
        type_labels.extend([atom_types[int(i-1)] for j in range(type_counts[int(i-1)])])
        
        
    configs: list[Atoms] = []
    for i in range(num_snapshots):
        config = Atoms(symbols=type_labels,positions=all_atomic_positions[i],cell=all_lattice_vectors[i],pbc=True)
        configs.append(config)
        
    return configs
        
        
        
        
        

class buildLAMMPS():
    '''
    build a set of LAMMPS calculation directories. 
    '''
    
    def __init__(self, base_directory: str, new_dir_label: str='lammps'):
        
        print('BASE:', base_directory)
        
        self.base_directory = base_directory
        
        self.lmps_input = list(Path(self.base_directory).glob('in.*'))
        self.lmps_dat = list(Path(self.base_directory).glob('*.dat'))
        
        assert self.lmps_input, 'Error: LAMMPS base directory must contain an input file formated as "in.*"'
        
        assert self.lmps_dat, 'Error: LAMMPS base directory must contain initial positions *.dat file'
        
        assert len(self.lmps_input) == 1, 'Error: base directory has multiple files named "in.* '
        
        self.var_in_file = open(self.lmps_input[0],'r').read()
        
        self.input_file_name = os.path.basename(self.lmps_input[0])
        
        # remove comments:
        var_in_no_com = []
        for line in self.var_in_file.splitlines():
            line_without_comment = line.split('#', 1)[0].rstrip()
            var_in_no_com.append(line_without_comment)
            
        self.var_in_file = "\n".join(var_in_no_com)
        
        lmps_input = list(Path(self.base_directory).glob('in.*'))
        
        # label defines how output dirs are named
        self.label = new_dir_label 
        
        pass
    
    def fetch_variables(self):
        '''
        Searches for all variables.
        '''
        
        self.variables = list(set(re.findall(r"£\S+", self.var_in_file)))
        
        self.variables_dict = {var: None for var in self.variables}
        
        for line in self.var_in_file.splitlines():
            
            linesplit = line.split()
            
            if linesplit:
                if linesplit[0] in self.variables:
                    
                    # relies on format: £variable = 
                    # since '=' will be second in the list
    
                    self.variables_dict[linesplit[0]] = linesplit[1:]
                    

        assert None not in self.variables_dict.values(), 'You called an undefined a variable'

    def set_variables(self):
        
        self.new_input_files = []
        
        self.new_dir_names = []
        
        all_values = list(self.variables_dict.values())
        
        # If multiple variables we need to specify all possible combinations 
        
        for combination in product(*all_values):
            
            new_dir_name = self.label
            
            current_variables_dict = {var: None for var in self.variables}
            
            for i, val in enumerate(combination):
                
                current_variables_dict[self.variables[i]] = val
                
                var_name = self.variables[i].replace('£','')
                new_dir_name += f'_{var_name}_{val}'
                
            new_input_file_str = []
            
            for line in self.var_in_file.splitlines():
            
                linesplit = line.split()
                
                pass_line = False
                
                if linesplit:
                    
                    for var in self.variables:
        
                        if linesplit[0] == var:
                            pass_line = True
                             # we want to ignore these lines for out final output as lammps will fail if line starts with £
                        
                        if var in linesplit:
                            
                            linesplit = [current_variables_dict[var] if x == var else x for x in linesplit]
                
                if pass_line:
                    continue 
                
                
                new_line = " ".join(linesplit)
                             
                new_input_file_str.append(new_line)
                
            new_input_file_str = "\n".join(new_input_file_str)
            
            self.new_input_files.append(new_input_file_str)
            
            self.new_dir_names.append(new_dir_name)
        
        
    def generate_outputs(self, output_directory: str='.'):
        '''
        in the output directory the code will generate: output_directory/dir1 output_directory/dir2 etc
        '''
        
        self.lammps_output_dir = output_directory
        
        Path(output_directory).mkdir(exist_ok=True)
        
        self.new_dirs = []
        
        for i, input_file in enumerate(self.new_input_files):
            
            new_dir = output_directory + '/' + self.new_dir_names[i]
            
            self.new_dirs.append(new_dir)
            
            shutil.copytree(self.base_directory, new_dir ,dirs_exist_ok=True)
            
            with open(output_directory + '/' + self.new_dir_names[i] + '/' + Path(self.lmps_input[0]).name, 'w') as f:
                
                f.write(input_file)
                
            print(f'Caculation directory generated: {output_directory}/{self.new_dir_names[i]}')
    
    def generate_bash_command(self,lammps_cmd_line: str):
        
        dirs = ''
        for dir in self.new_dir_names:
            dirs += f'{dir} '
            
        output_name = 'out.' + self.input_file_name.split(".")[1]
        
        
        result = f'directories="{dirs}"\n'
        
        result += f'cd {self.lammps_output_dir}\n'
        
        result+=f'''for i in $directories; 
do 
cd $i
    {lammps_cmd_line} -i {self.input_file_name} -l {output_name}
cd ..
done\n'''
        result += 'cd ..'

        
        return result
        
            
    

    
    
    

#base_dir = sys.argv[1]
#out_dir = sys.argv[2]

#test = contructLAMMPS(f'./{base_dir}')
#test.fetch_variables()
#test.set_variables()
#test.generate_outputs(output_directory=f"./{out_dir}")