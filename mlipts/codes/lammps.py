'''
@author: William Davie

File containing lammps specific functionality. Used to build many lammps calculations.

Main functionality is to automatically set up many LAMMPS calculations, snapshots of each calculation can then used for QM caclulations. 

The user provides a 'base directory' including:

    The input (in.*) file:
        - this should include the potential specification, often this will refer to a file elsewhere and a minimization task. 
        - Variables marked with '£' which will be replaced by values in a list/array.
        
    The .dat file containing atomic positions
'''

import os, sys
import shutil
from pathlib import Path
import numpy as np
import re
from itertools import product
from ase import Atoms
from typing import Optional 


def read_lammps_output(outdir: str, atom_types: list[str], pbc: Optional[bool]=True) -> list[Atoms]:
    '''
    Read lammps *md output file.
    
    Parameters
    ----------
    outdir: str
        directory of lammps calculation
    atom_types: str[list]
        list of atom types, ordered as in lammps output
    pbc: Optional[bool]
        define the type of periodic boundry conditions used in simulation.
    
    Returns 
    -------
    configs: list[Atoms]
        list of atomic configurations fetched from dump file. Returned as a list of ase.Atoms types.
    '''
    
    outdir = list(Path(outdir).glob('md.*'))
    md_output = open(outdir[0],'r').read()
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


def write_lammps_submission_script(lammps_dirs: list[str], lammps_cmd_line: str, output_directory: str ='.'):
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
cd -
done\n'''
    result += 'cd ..'
    
    return result


def build_lammps_calculations(base_dir: str, variables: dict, outdir: str='.') -> list[str]:
    '''
    Generates a set of directories containing input files for lammps calculations. 
    Reads the base directory, searches for a set of marked variables '£' and generates a directory for all values in the provided dictionary.
    
    Parameters
    ----------
    base_dir: str
        path to the base lammps directory containing a in.* and .dat file. in.* is expected to include some variables marked '£'
    variables: dict
        keys should be formated as '£*' with corrosponding arrays.
    outdir: str 
        output directory for new calculations  
        
    Returns
    -------
    new_dirs: list[str]
        a list of paths containing all calculations to be run.
        
    Raises
    ------
    FileNotFoundError:
        if in.* or *.dat not found in base_dir
    NameError:
        if a key in variables cannot be found in *.in
    
    '''
    
    build = lammpsBuild(f'{base_dir}',variables)
    build.read_base_directory()
    build.generate_calculations(outdir=outdir)
    return build.new_dirs


class lammpsBuild():
    
    def __init__(self, lammps_base_dir: str, variables: dict) -> None:
        '''
        Initialize a lammps build. 
        
        Parameters
        ----------
        lammps_base_dir: str
            path to directory containing files for a lammps calculation
        variables: dict
            dictionary of variables to build new lammps directories from.
        '''
        self.lammps_base_dir = lammps_base_dir
        
        self.input = None
        
        self.variables = variables
        self.variable_keys = list(self.variables.keys())
        
        self.new_dirs = []
        
        return None

        
    def read_base_directory(self) -> None:
        '''
        Reads files and errors if base directory does not have the correct format for constructing multiple directories.
        
        Raises
        ------
        FileNotFoundError:
            if in.* or *.dat not found in base_dir
        FileExistsError:
            if multiple in.* or *.dat files are found
        NameError:
            if a key in variables cannot be found in *.in
        '''
        
        self.input_files = list(Path(self.lammps_base_dir).glob('in.*'))
        self.dat_files = list(Path(self.lammps_base_dir).glob('*.dat'))
        
        if not self.input_files:
            raise FileNotFoundError('No input (in.*) file found in base directory')
        if not self.dat_files:
            raise FileNotFoundError('No atomic data (*.dat) file found in base directory')
        if len(self.input_files) > 1:
            raise FileExistsError('Cannot have multiple input files (in.*)')
        if len(self.dat_files) > 1:
            raise FileExistsError('Cannot have multiple atomic data files (*.dat)')

        self.input = open(self.input_files[0],'r').read()
    
        # remove comments
        in_no_com = []
        for line in self.input.splitlines():
            line_without_comment = line.split('#', 1)[0].rstrip()
            in_no_com.append(line_without_comment)
        self.input = "\n".join(in_no_com)
        
        variables_in_file = list(set(re.findall(r"£\S+", self.input)))

        for i in self.variable_keys:
            if i not in variables_in_file:
                raise NameError(f'Variable {i} in input dictionary not found in lammps input file (in.*)')
        
        return None
    
    
    def generate_calculations(self, label: str='lammps', outdir: str='.') -> None:
        '''
        Given variables and a base directory, a set of new calculation directories are generated.
        '''
        
        new_input_files = []
        new_dir_names = []
        all_values = list(self.variables.values())
        
        if self.input is None:
            self.read_base_directory()
        
        # If multiple variables we need to specify all possible combinations 
        for combination in product(*all_values):
            #combination looks like tuple(var1, var2, var2, etc) 
            new_dir_name = label
            new_input_str = []
            current_variables = {var: None for var in self.variable_keys}
            
            for i, val in enumerate(combination):
                # label new calculation directory
                var_name = self.variable_keys[i][1:] # assumes first character is £
                new_dir_name += f'_{str(var_name).replace("£","")}_{round(val,None)}'
                current_variables[self.variable_keys[i]] = val
                
            #now define new str for this combination 
            for line in self.input.splitlines():
                line_split = line.split()
                
                for var_name in self.variable_keys:
                    if var_name in line_split:
                        line_split = [str(current_variables[var_name]) if x == var_name else x for x in line_split]
                        
            
                new_line = " ".join(line_split)
                new_input_str.append(new_line)
                
            new_input_str = "\n".join(new_input_str)
            
            self.__write_calculation__(new_dir_name,new_input_str,outdir)
            
        return None
            
                
    def __write_calculation__(self, calc_name: str, input_str: str, outdir: str) -> None:
        '''
        Given a base directory (self), calculation name and lammps input file string, a new calculation is written.
        '''
        
        Path(outdir).mkdir(exist_ok=True)
    
        new_dir = outdir + '/' + calc_name 
        self.new_dirs.append(new_dir)
        shutil.copytree(self.lammps_base_dir, new_dir ,dirs_exist_ok=True)
            
        with open(new_dir + '/' + Path(self.input_files[0]).name, 'w') as f:
            f.write(input_str)
                
        print(f'Caculation directory generated: {outdir}/{new_dir}')
        
        return None