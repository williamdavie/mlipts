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
        
            
        
    
  
  
def build_lammps_calculations(base_dir: str,lammps_cmd_line: str=None) -> list[str]:
    '''
    Generates a set of directories containing input files for lammps calculations
    
    Number of directories generated depends on the variable set in ./lammps_base_*

    Returns a list of directory names.
    '''
    
    build = buildLAMMPS(f'{base_dir}')
    build.fetch_variables() 
    build.set_variables()
    build.generate_outputs()
    if lammps_cmd_line == None:
        pass
    else:
        build.generate_bash_script(lammps_cmd_line=lammps_cmd_line)
    
    

    return build.new_dirs
    
    
if __name__ == '__main__':
    
    test = buildLAMMPS(f'./base_example')
    test.fetch_variables()
    test.set_variables()
    test.generate_outputs('./Examples')
    test.generate_bash_script('lmps')
    
    
    

#base_dir = sys.argv[1]
#out_dir = sys.argv[2]

#test = contructLAMMPS(f'./{base_dir}')
#test.fetch_variables()
#test.set_variables()
#test.generate_outputs(output_directory=f"./{out_dir}")



