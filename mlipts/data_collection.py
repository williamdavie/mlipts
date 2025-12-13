
'''
@author William Davie

Main class for collecting a data set for model training/fine-tuning. 
'''

from mlipts.codes.lammps import build_lammps_calculations, read_lammps_output
from mlipts.codes.vasp import build_vasp_calculation
from mlipts.hpc_submission.archer2 import archer2_submission_template
from mlipts.similarity.filter import filter_by_emd
from ase import Atoms
from ase.io import read, write
from pathlib import Path
import subprocess
from mlipts.constants import __hpcs__,__diffmethods__,__MDcodes__,__QMcodes__


class DataCollection():
    
    def __init__(self, atom_types: list[str]):
        '''
        Initialize a DataCollection class. Used to develop a database for training a Machine Learned Interatomic Potential.
        '''

        self.atom_types = atom_types
        
        '''
        Both MD and QM calculations are stored in states: initialized, active, complete
        
        initialized: directory created but not run
        active: simulation complete, active for usage.
        complete: simulation complete, inactive.
        '''
        

        self.initialized_MD_dirs = []
        self.active_MD_dirs = [] # simulation in directory run, active for analysis
        self.active_MD_configs: list[Atoms] = []
        self.complete_MD = [] # simulation stored as complete, no longer active
        
        self.initialized_QM_dirs = []
        self.active_QM_dirs = []
    
        
        self.MD_submission_count = 0

        
    def build_MD_calculations(self, 
                              MD_base_dir: str, 
                              variables: dict,
                              MDcode: str='lammps',
                              outdir: str='.') -> None:
        '''
        Generates a set of directories for Molecular Dynamics simulations given the parameters set in MD_base. 
        
        Parameters
        ----------
        MD_base_dir : str
            A directory containing all necessary files for successfully running a simulation given an MD code. See docs for details.
        MDcode: str 
            MD code of choice. Default is lammps.
        outdir: str
            Define a directory to store calculation files. Default is the working directory.
            
        Returns
        -------
        None : None
            New directories generated in working directory. Directory paths added to self.MD_calculations
            
        Raises
        ------
        ValueError
            if chosen MD code not supported. 
        '''
        
        if MDcode == 'lammps':
            new_dirs = build_lammps_calculations(MD_base_dir,variables,outdir=outdir)
            self.initialized_MD_dirs.extend(new_dirs)
        elif MDcode not in __MDcodes__:
            raise ValueError(f'{MDcode} not supported.')
        
        return None
    
    def write_MD_submission_scripts(self, MD_cmd_line: str,
                                   nodes: int, ranks: int,
                                   time: str,
                                   MDcode: str='lammps',
                                   hpc: str='archer2',
                                   npartitions: int=1,
                                   scripts_outdir: str='./MD_scripts',
                                   submit: bool=True,
                                   mark_as_active: bool=True,
                                   header_str: str=None,
                                   hpc_account: str=None):
        '''
        write submission script for Molecular Dynamics simulations, built for all directories marked 'initialized'. 
        
        Parameters
        ----------
        MD_cmd_line: str
            the command line used to run the chosen MD code (see examples).
        nodes: int
        ranks: int
        time: str
            run time formated as "XX:XX:XX"
        MDcode : str
            MD code of choice
        hpc : str
            hpc of choice for header of submission script.
        submit: bool
            automatically sbatch the script. 
        
        Returns
        -------
        None : None
            generates MD_submission_script_#i in the working directory.
        
        '''
        # header
        header = fetch_hpc_header(hpc,header_str,nodes,ranks,time,hpc_account)
        
        # cmds
        if MDcode not in __MDcodes__:
            raise ValueError(f'MD code {MDcode} not supported.')
        
        cmd_scipts = write_run_calculation_scripts(self.initialized_MD_dirs,MD_cmd_line,npartitions=npartitions) # leaving as one partition as default for now but more partitions possible easy addition.
        Path(scripts_outdir).mkdir(exist_ok=True)
        for i,cmd in enumerate(cmd_scipts):
            with open (f'{scripts_outdir}/MD_submission_script_#{i}','w') as f:
                f.write(header)
                f.write('\n')
                f.write(cmd)
                
            print(f'MD submission script saved to: {scripts_outdir}/MD_submission_script_#{i}')
            
            if submit:
                subprocess.run(f'sbatch {scripts_outdir}/MD_submission_script_#{i}',shell=True)
        
        if mark_as_active:
            self.active_MD_dirs.extend(self.initialized_MD_dirs)
            self.initialized_MD_dirs = []

      
        return None
    

    def filter_active_MD(self,
                         tol: float, 
                         method: str='emd',
                         auto: bool=False, 
                         k: int=20, 
                         show_dendrograms: bool=False) -> None:
        '''
        Removes some configurations from set of active MD configs if they are too similar to eachother. The similarity between configurations is caclulated using a distance metric defined by method.
        
        Parameters
        ----------
        tol: float 
            if the difference metric between two configurations is lower than tol, one of the configurations will be removed.
        method: str  
            method to calculate a difference metric. Default is earth movers distance "emd".
        auto: bool
            auto select command line input options. Default is False.
            
        Raises
        ------
        ValueError:
            if difference metric method not availible.
        '''
        
        if not self.active_MD_configs:
            self.fetch_MD_configs_from_calcs()
        
        if method == 'emd':
            new_configs, inds = filter_by_emd(self.active_MD_configs,tol,k=k,show_dendrograms=show_dendrograms)
        elif method not in __diffmethods__:
            raise ValueError(f'Distance metric method {method} not found')
    
        print(f"Filter reduced configuration space from {len(self.active_MD_configs)} to {len(new_configs)}")
        
        update_active = input("Update active configs? ('Y'/'N'): ")
        if update_active.capitalize() == 'Y':
            self.active_MD_configs = []
            self.active_MD_configs.extend(new_configs)
            print("Active configs updated.")
    
    
        return None
    
    def build_QM_calculations(self, 
                              QM_base_dir: str, 
                              QMcode: str='vasp', 
                              outdir: str = './QM_calculations', 
                              label: str='') -> None:
        '''
        For all active MD configs (self.active_MD_configs) a first principle calculation directory is generated.
        
        Parameters
        ----------
        QM_base_dir: str
            directory containing necessary files for a QM calculations except atomic position information, which will be taken from active MD calculations.
        QMcode: str
            code used. Default is vasp.    
            
        Returns
        -------
        None: None
            new QM calculations generated in outdir. 
        '''
        
        if not self.active_MD_configs:
            self.fetch_MD_configs_from_calcs()
            
        print(f'Number of active configs = number of QM calculation directories = {len(self.active_MD_configs)}')
        
        for i,config in enumerate(self.active_MD_configs):
            if QMcode == 'vasp':
                new_calc_dir = build_vasp_calculation(QM_base_dir,config,f'{label}_c_#{i}',outdir)
            elif QMcode not in __QMcodes__:
                raise ValueError(f'QM code {QMcode} not supported')
            
            self.initialized_QM_dirs.append(new_calc_dir)
            
        print(f'Calculations stored in {outdir}')
        
        return None
    
    def write_QM_submission_scripts(self, QM_cmd_line: str,
                                   nodes: int, ranks: int,
                                   time: str,
                                   npartitions: int=1,
                                   save_and_remove: bool=True,
                                   QMcode: str='vasp',
                                   python_env: str=None,
                                   database_file: str=None,
                                   hpc: str='archer2',
                                   hpc_account: str=None,
                                   scripts_outdir: str='./QMscripts',
                                   submit: bool=True,
                                   header_str: str=None,
                                   mark_as_active: bool=True):
        '''
        Build submission script for Quantum Mechanical (first principle) simulations, built for all directories marked 'initialized'. 
        
        Parameters
        ----------
        QM_cmd_line: str
            the command line used to run the chosen QM code (see examples).
        nodes: int
        ranks: int
        time: str
            run time formated as "XX:XX:XX"
        QMcode : str
            QM code of choice
        hpc : str
            hpc of choice for header of submission script.
        python_env: str
            path to python enviroment
        database_file: str
            path to file to store data.
        
        
        Returns
        -------
        None : None
            generates QM_submission_script_#i in the working directory.
            
        Raises
        ------
        see write_run_calculation_scripts()
        
        '''
        
        # header
        header = fetch_hpc_header(hpc,header_str,nodes,ranks,time,hpc_account)
        
        #cmd
        if QMcode not in __QMcodes__:
            raise ValueError(f'QM code {QMcode} not supported')
        cmd_scipts = write_run_calculation_scripts(self.initialized_QM_dirs,
                                                   QM_cmd_line,
                                                   npartitions=npartitions,
                                                   save_and_remove=save_and_remove,
                                                   python_env=python_env,code=QMcode,database_file=database_file)
        Path(scripts_outdir).mkdir(exist_ok=True)
        for i,cmd in enumerate(cmd_scipts):
            with open (f'{scripts_outdir}/QM_submission_script_#{i}','w') as f:
                f.write(header)
                f.write('\n')
                f.write(cmd)
            print(f'MD submission script saved to: {scripts_outdir}/QM_submission_script_#{i}')
            if submit:
                subprocess.run(f'sbatch {scripts_outdir}/QM_submission_script_#{i}',shell=True)
        
        
        if mark_as_active:
            self.active_MD_dirs.extend(self.initialized_MD_dirs)
            self.initialized_MD_dirs = []
    
        return None

    
    def fetch_MD_configs_from_calcs(self):
        '''
        Collect the atomic configs from the active MD directories
        '''
        
        if not self.active_MD_dirs:
            raise ValueError('Tried to collect configurations from active MD directories but no directories are active.')
        
        for dir in self.active_MD_dirs:
            self.active_MD_configs.extend(read_lammps_output(dir,atom_types=self.atom_types))
        return self.active_MD_configs
    
    
    def set_active_MD_dirs(self, outdir: str='./MD_calculations') -> None:
        '''
        Set MD calculation directories manually.
        '''
        
        old_len = len(self.active_MD_dirs)
        path = Path(outdir)
        subdirs = [p for p in path.iterdir() if p.is_dir()]
        for calc in subdirs:
            if calc not in self.active_MD_dirs:
                self.active_MD_dirs.append(str(calc))
                
        print(f'Number of active MD directories updated from {old_len} to {len(self.active_MD_dirs)}')
        
        return None
    
    def set_active_MD_configs(self,
                              config_file: str) -> None:
        '''
        Adds configurations in config_file to self.active_MD_configs.
        '''
        
        old_len = len(self.active_MD_configs)
        configs = read(config_file,':')
        self.active_MD_configs.extend(configs)
        print(f'Number of active configs updated from {old_len} to {len(self.active_MD_configs)}')
        
        return None

    
    def __check_initialized__(self):
        '''
        prints initialized MD directories
        '''
        print('MD directories initialized: ')
        print(self.initialized_MD_dirs)
        
    
    def check_active_MD(self):
        '''
        prints active MD directories
        '''
        print('MD directories currently active: ')
        print(self.active_MD_dirs)
        print('Num active: ', len(self.active_MD_dirs))
        
    def __save_active__(self, outname: str):
        
        write(f'{outname}.xyz',self.active_MD_configs)
        
        
def fetch_hpc_header(hpc: str, header_str: str, 
                             nodes: int=1,
                             ranks: int=1,
                             time: str='01:00:00',
                             hpc_account: str=None) -> str:
    '''
    Given some hpc parameters returns the header of a submission script. 
    '''

    if hpc == 'archer2':
        header = archer2_submission_template(nodes,ranks,time,account=hpc_account)
    elif hpc == 'custom':
        if header_str == None:
            raise ValueError('custom hpc header but no header_str argument provided.')
        else:
            header = header_str
    elif hpc not in __hpcs__:
        raise ValueError(f'hpc {hpc} not supported.')
    
    return header
        
        
        
def write_run_calculation_scripts(calc_dirs: list[str],
                                 cmd_line: str,
                                 npartitions: int=1, 
                                 save_and_remove: bool=False,
                                 python_env: str=None,
                                 code: str=None,
                                 database_file: str=None) -> list[str]:
    '''
    Given a list of calculation directories, npartitions scipt(s) are generated to enter each directory and run a command line.
    
    Parameters
    ----------
    calc_dirs: list[str]
        list of directories to enter and run the cmd_line
    cmd_line: str
        command line to run in each directory
    npartitions : int
        number of scripts to generate. Calculations per script is len(calc_dirs/npartitions)
    save_and_remove: bool
        save data into a configuration database and remove the calculation outputs, used to save disk space.
    savedata_cmd: str
        command run to save data from a given calculation. 
    
    Returns
    -------
    calc_scripts : list[str]
        list of paths to each submission script
        
    Raises
    ------
    ValueError
        if number of calculations is not divisible by number of paritions. 
    ValueError
        if save_and_remove is set to true but no command is given.
    FileNotFoundError
        if python_env doesn't include /bin/python
    '''
    
    calc_scripts = []
    
    if len(calc_dirs) % npartitions != 0:
        raise ValueError(f'Number of calculations to run ({len(calc_dirs)}) is not divisible by number of partitions specified {npartitions}')

    if save_and_remove == True:
        if python_env == None or code == None or database_file == None:
            raise ValueError(f'save_and_remove option requires specification of the code used, database to save to and a python enviroment with mlipts.')
        if not Path(f'{python_env}/bin/python').exists():
            raise FileNotFoundError(f"couldn't find python at ({python_env}/bin/python)")
        
        savedata_cmd = f'{python_env}/bin/python -m mlipts.append_to_database $i {database_file} {code}'
        remove_cmd = 'rm -r $i'
    else:
        savedata_cmd = ''
        remove_cmd = ''
    
    num_calcs_per_submission = int(len(calc_dirs) / npartitions)

    for i in range(npartitions):
        
        script = ''

        current_dirs=''
        for dir in calc_dirs[int(i*num_calcs_per_submission):int((i+1)*num_calcs_per_submission)]:
            current_dirs+=f'{dir} '
        
        script+=f'directories="{current_dirs}"\n'
        script+=f'''for i in $directories; 
do 
cd $i
{cmd_line}
cd -
{savedata_cmd}
{remove_cmd}
done\n'''

        calc_scripts.append(script)
            
    return calc_scripts
            