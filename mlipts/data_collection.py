
'''
@Author William Davie

Main class for collecting a data set for model training/fine-tuning. 
'''

from mlipts.codes.lammps import *
from mlipts.codes.vasp import *

from mlipts.hpc_submission.archer2 import *

from mlipts.append_to_database import *

class trainingWorkflow():
    
    def __init__(self, MD_base: str, electronic_base: str):
        
        self.MD_base = MD_base
        self.electronic_base = electronic_base
        
        # stores the script that runs MD calculations for our workflow. 
        self.MD_sequence_script = None
        
        # stores a list of directories containing MD calculations
        self.MD_calculations = None
        self.electronic_calculations = None
        
        # a list of directory strings
        self.electronic_calculation_dirs: list[str] = []
        

    def build_lammps_calculations(self) -> list[str]:
        '''
        Generates a set of directories containing input files for lammps calculations
    
        Number of directories generated depends on the variable set in ./lammps_base_*

        Returns a list of directory names.
        '''
        
        build = buildLAMMPS(f'{self.MD_base}')
        build.fetch_variables() 
        build.set_variables()
        build.generate_outputs()

        self.MD_calculations = build
        
    def build_lammps_submission_script(self,lammps_cmd_line: str,nodes: int, ranks: int, time: str, hpc: str='Archer2', show: bool=False):
        '''
        Writes a script to run 
        '''
        
        assert type(self.MD_calculations) == buildLAMMPS, 'Your MD calculation is not LAMMPS format, must run build_lammps_calculations()'
        
        lammps_script = ''
        
        
        if hpc == 'Archer2':
            
            lammps_script += archer2_submission_template(nodes=nodes,ranks=ranks,time=time)
        
        lammps_script+='\n'
        
        lammps_script += self.MD_calculations.generate_bash_command(lammps_cmd_line)
        
        
        
        if show==True:
            
            print('SCRIPT PREVIEW: ')
            print(lammps_script)
            
        with open('submit_lammps_calcs','w') as f:
            f.write(lammps_script)
            print('Submission script saved to: submit_lammps_calcs')
            
    
    def build_vasp_for_lammps(self,atom_type_labels: list[str]):
        
        assert type(self.MD_calculations) == buildLAMMPS, 'Your MD calculation is not LAMMPS format, must run build_lammps_calculations()'
        
        self.electronic_calculations: list[constructVASP] = []
        # for all lammps calculations in current build
        self.DFTcount: int = 0
        
        for lammps_calc in self.MD_calculations.new_dirs:
        
            build_vasp = constructVASP(self.electronic_base,lammps_calc)
            build_vasp.read_pos_LAMMPS()
            build_vasp.build_vasp_dirs(atom_type_labels)
            
            self.electronic_calculations.append(build_vasp)
            
            for dir in build_vasp.new_vasp_calc_dirs:
                
                self.electronic_calculation_dirs.append(dir)
            
            self.DFTcount += len(build_vasp.new_vasp_calc_dirs)
        
        # this is currently not quite correct. 
        print(f' DFT calculations to run with {build_vasp.num_atoms} atoms: {self.DFTcount}')
        
    def build_vasp_submission_scripts(self, vasp_cmd_line: str, num_script_partitions: int, nodes: int, ranks: int, time_per_script: str, hpc: str='Archer2', show: bool=False):
        
        # Leaving this as is for now:
        assert self.DFTcount % num_script_partitions == 0, f'Number DFT counts ({self.DFTcount}) must be divisible by number of partions ({num_script_partitions})'
        num_calcs_per_submission = self.DFTcount / num_script_partitions

        '''
        For DFT calculations it is a lot faster to partition calculations over a number of submission scripts. This is done via: num_script_partitions.
        '''
        
        print(f'Each script will carry out {num_calcs_per_submission} DFT caculations, over a time of {time_per_script}')
        
        
        for i in range(num_script_partitions):
            
            vasp_submit = ''
        
            if hpc == 'Archer2':
            
                vasp_submit += archer2_submission_template(nodes=nodes,ranks=ranks,time=time_per_script)
        
            vasp_submit+='\n'
            
            current_dirs = ''
            for dir in self.electronic_calculation_dirs[int(i*num_calcs_per_submission):int((i+1)*num_calcs_per_submission)]:
                current_dirs += f'{dir} '
            
            vasp_submit += f'directories="{current_dirs}"\n'

            vasp_submit+=f'''for i in $directories; 
do 
cd $i
    {vasp_cmd_line}
cd -
done\n'''
            
            if show==True:
                print('SCRIPT PREVIEW: ')
                print(vasp_submit)
                
            with open(f'submit_vasp_#{i}','w') as f:
                f.write(vasp_submit)
                print(f'Submission script saved to: submit_vasp_#{i}')

        
    
    def append_vasp_to_mace_data(self,database_file: str, all: bool=True, clean_vasp_dirs: bool=True):
        '''
        For all the vasp calculations in this workflow, collect the data. 
        
        The data will append to the input dataset_file
        '''
        
        for dft_calc in self.electronic_calculation_dirs:
            
            # periodic boundry conditions are default here 
            append_vasp_calculation(dft_calc,database_file)
            
    
        return None
        

    #def build_archer2_MD_submission(self):
        
        

    