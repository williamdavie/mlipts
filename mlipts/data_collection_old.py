
'''
@author William Davie

Main class for collecting a data set for model training/fine-tuning. 
'''

from mlipts.codes.lammps import *
from mlipts.codes.vasp import *
from mlipts.hpc_submission.archer2 import *
from mlipts.append_to_database import *
import subprocess

class DataCollection():
    
    def __init__(self, MD_base: str, electronic_base: str, atom_types: list[str]) -> None:
        
        self.MD_base = MD_base
        self.electronic_base = electronic_base
        
        # stores the script that runs MD calculations for our workflow. 
        self.MD_sequence_script = None
        
        self.QM_scripts: list[str] = []  # a list of submission scripts for DFT.
        
        # stores a list of directories containing MD calculations
        self.MD_calculations = None
        self.electronic_calculations = None
        
        # a list of directory strings
        self.electronic_calculation_dirs: list[str] = []
        
        self.atom_types = atom_types
        

    def build_lammps_calculations(self,label: str='lammps') -> list[str]:
        '''
        Generates a set of directories containing input files for lammps calculations
    
        Number of directories generated depends on the variable set in MD_base

        Returns a list of directory names.
        '''
        
        build = buildLAMMPS(f'{self.MD_base}',label)
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
            
    
    def build_vasp_for_lammps(self):
        
        assert type(self.MD_calculations) == buildLAMMPS, 'Your MD calculation is not LAMMPS format, must run build_lammps_calculations()'
        
        self.electronic_calculations: list[constructVASP] = []
        # for all lammps calculations in current build
        self.DFTcount: int = 0
        
        for lammps_calc in self.MD_calculations.new_dirs:
            
            md_output = list(Path(lammps_calc).glob('md.*'))

            if not md_output:
        
                print(f'''Error: Lammps simulation directory ({lammps_calc}) does not contain output file formated as "md.*"
                
        - This calculation will not be used to construct VASP directories.
        - To use this LAMMPS calculation, rerun this method. 
                ''')
    
                continue
    
        
            build_vasp = constructVASP(self.electronic_base,lammps_calc)
            build_vasp.read_pos_LAMMPS()
            build_vasp.build_vasp_dirs(self.atom_types)
            
            self.electronic_calculations.append(build_vasp)
            
            for dir in build_vasp.new_vasp_calc_dirs:
                
                self.electronic_calculation_dirs.append(dir)
            
            self.DFTcount += len(build_vasp.new_vasp_calc_dirs)
        
        # this is currently not quite correct. 
        print(f' DFT calculations to run with {build_vasp.num_atoms} atoms: {self.DFTcount}')
        
    def build_vasp_submission_scripts(self, vasp_cmd_line: str, num_script_partitions: int, nodes: int, ranks: int, time_per_script: str, hpc: str='Archer2', show: bool=False, 
                                      save_and_remove: bool=False, database_file: str=None, pythonenv: str=None):
        
        # Leaving this as is for now
        assert self.DFTcount % num_script_partitions == 0, f'Number DFT counts ({self.DFTcount}) must be divisible by number of partions ({num_script_partitions})'
        num_calcs_per_submission = self.DFTcount / num_script_partitions

        '''
        For DFT calculations it is a lot faster to partition calculations over a number of submission scripts. This is done via: num_script_partitions.
        '''
        
        print(f'Each script will carry out {num_calcs_per_submission} DFT caculations, over a time of {time_per_script}')
        
        
        if save_and_remove == True:
            assert database_file != None, "Error you have save data as your calculations are performed but did not specify a database_file name"
            assert pythonenv != None, 'Error you have save data, this requires python, please define a python enviroment.'
            # relies on a lot of formatting across the code
            savedata_cmd = f'{pythonenv}/bin/python -m mlipts.append_to_database $i {database_file}'
            remove_cmd = 'rm -r $i'
            
        else:
            savedata_cmd = ''
            remove_cmd = ''
        
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
{savedata_cmd}
{remove_cmd}
done\n'''
            
            if show==True:
                print('SCRIPT PREVIEW: ')
                print(vasp_submit)
            
            Path('./QM_submission').mkdir(exist_ok=True)
            
            with open(f'./QM_submission/submit_vasp_#{i}','w') as f:
                f.write(vasp_submit)
                self.QM_scripts.append(f'./QM_submission/submit_vasp_#{i}')
                
                print(f'Submission script saved to: ./QM_submission/submit_vasp_#{i}')

    
    def QM_submit_all(self, start: int=0, end: int=None) -> None:
        '''
        sbatch ~all scripts~
        
        start and end define the index of your scripts. Some hpcs will limit how many you can submit.
        '''
        
        print('This command must be run from your working directory.')
        
        # save IDs if need to cancel job. 
        self.jobIDs = []
        
        if end==None: end=len(self.QM_scripts)
        
        for i in self.QM_scripts[start:end]:
            
            result = subprocess.run(f'sbatch {i}', shell=True, capture_output=True, text=True)
            
            result = result.stdout
            print(result)
            self.jobIDs.append(result.split()[2])
            
    def QM_cancel_all(self, start: int=0, end: int=None):
        '''
        cancel all scripts if you made a mistake using QM_submit_all()
        '''
        if end==None: end=len(self.QM_scripts)
    
        for i in self.QM_scripts[start:end]:
            
            print(i)
            
            result = subprocess.run(f'scancel {i}', shell=True, capture_output=True, text=True)
            
            print(result.stdout)
            

    def append_lammps_to_dataset(self, dataset: str):
        '''
        Adds lammps configurations to a dataset, used for similarity asse
        '''
        
        return None
        
    
    def append_vasp_to_mace_data(self,database_file: str, all: bool=True, clean_vasp_dirs: bool=True):
        '''
        For all the vasp calculations in this workflow, collect the data. 
        
        The data will append to the input dataset_file
        '''
        
        for dft_calc in self.electronic_calculation_dirs:
            
            # periodic boundry conditions are default here 
            append_vasp_calculation(dft_calc,database_file)
            
    
        return None
    
    def clean_all(self):
        '''
        Deletes all calculation directories.
        '''
        
        print('Warning <!>, you are about to delete all calculation folders, ensure your data has been saved')
        continue_bool = input('Are you sure you want to continue "Y" / "N": ')
        
        if continue_bool == 'Y':
            
            for dir in self.MD_calculations.new_dirs:
                try:
                    shutil.rmtree(dir)
                except:
                    pass
                



    #def build_archer2_MD_submission(self):
        
        

    
def read_pos_LAMMPS(lammps_dir: str):
    '''
    From file md.* reads output as specified in LAMMPS. 
    '''
    
    md_output = open(lammps_dir).read()
    split_lines = md_output.splitlines()
    
    timesteps = [] # snapshot time steps.
    all_atomic_positions = [] 
    all_lattice_vectors = []
    num_snapshots = md_output.count('ITEM: TIMESTEP')
    num_atoms = int(split_lines[split_lines.index('ITEM: NUMBER OF ATOMS') + 1])
    
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