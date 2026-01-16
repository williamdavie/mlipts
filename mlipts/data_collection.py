
'''
@author William Davie

Main class for collecting a data set for model training/fine-tuning. 
'''

import numpy as np 
import ase
import ase.io
from pathlib import Path
import subprocess
import yaml
from mlipts.codes import lammps, vasp
from mlipts.hpc_submission import hpc_utils
from mlipts.similarity import filter,group
from mlipts.constants import __hpcs__,__diffmethods__,__MDcodes__,__QMcodes__


class DataCollection():
    
    def __init__(self, atom_types: list[str],
                 hpc_config: str=None, **kwargs):
        '''
        Initialize a DataCollection class. Used to develop a database for training a Machine Learned Interatomic Potential.
        '''
        # hpc configuration for writing scripts.
        self.hpc_config = hpc_utils.load_hpc_config(self,hpc_config, **kwargs)
            
        '''
        Both MD and QM calculations are stored in states: initialized, active, complete
        
        initialized: directory created but not run
        active: simulation complete, active for usage.
        complete: simulation complete, inactive.
        '''

        self.initialized_MD_dirs: list[str] = []
        self.active_MD_dirs: list[str] = [] # simulation in directory run, active for analysis
        self.active_MD_configs: list[ase.Atoms] = []

        self.initialized_QM_dirs: list[str] = []
        self.active_QM_dirs: list[str] = []
        
        self.QM_base_dir = None
        
        self.submitted_jobs = []

        
    def build_MD_calculations(self, 
                              MD_base_dir: str, 
                              variables: dict,
                              MDcode: str='lammps',
                              outdir: str='MD_calculations') -> None:
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
            new_dirs = lammps.build_lammps_calculations(MD_base_dir,variables,outdir=outdir)
            self.initialized_MD_dirs.extend(new_dirs)
        elif MDcode not in __MDcodes__:
            raise ValueError(f'{MDcode} not supported.')
        
        return None
    
    def write_MD_submission_scripts(self, MD_cmd_line: str, time_per_partition: str,
                                   MDcode: str='lammps',
                                   npartitions: int=1,
                                   scripts_outdir: str='scripts',
                                   submit: bool=True,
                                   mark_as_active: bool=True,
                                   save_and_remove: bool=True,
                                   database_file: str='MD_samples.xyz',
                                   dependencies: list[str]=[]):
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
        
        self.hpc_config['time'] = time_per_partition # custom time.
        self.hpc_config['dependencies'] = dependencies
        header = hpc_utils.fetch_hpc_header(self.hpc_config)
        
        # cmds
        if MDcode not in __MDcodes__:
            raise ValueError(f'MD code {MDcode} not supported.')
        
        cmd_scipts = write_run_calculation_scripts(self.initialized_MD_dirs,
                                                   MD_cmd_line,
                                                   npartitions=npartitions,
                                                   save_and_remove=save_and_remove,
                                                   database_file=database_file,
                                                   code=MDcode,python_env=self.hpc_config['python_env']) # leaving as one partition as default for now but more partitions possible easy addition.
        Path(scripts_outdir).mkdir(exist_ok=True)
        for i,cmd in enumerate(cmd_scipts):
            script = hpc_utils.ScriptBuilder(header=header)
            script.add_cmd_line(cmd)
            script.write_script(f'{scripts_outdir}/MD_submission_script_#{i}',message=False)
            
            if submit:
                id = script.submit_script()
                self.submitted_jobs.append(id)
        
        if mark_as_active:
            self.active_MD_dirs.extend(self.initialized_MD_dirs)
            self.initialized_MD_dirs = []

      
        return None
    

    def filter_active_MD(self,
                         tol: float, 
                         method: str='emd',
                         auto: bool=False, 
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
            k = input('You opted to filter by earth movers distance, enter number of neighbours (k) to be considered: ')
            try: k = int(k)
            except: raise ValueError('k must be an interger.')
            new_configs, inds = filter.filter_by_emd(self.active_MD_configs,tol,k=k,show_dendrograms=show_dendrograms)
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
                              label: str=None,
                              pre_define_configs: list[ase.Atoms]=None) -> None:
        '''
        For all active MD configs (self.active_MD_configs) a first principle calculation directory is generated.
        
        Parameters
        ----------
        QM_base_dir: str
            directory containing necessary files for a QM calculations except atomic position information, which will be taken from active MD calculations.
        QMcode: str
            code used. Default is vasp. 
        pre_defined_configs
            
        Returns
        -------
        None: None
            new QM calculations generated in outdir. 
        '''
        
        if label==None:
            label=QMcode
        
        if not self.active_MD_configs and pre_define_configs is None:
            self.fetch_MD_configs_from_calcs()
            if not self.active_MD_configs:
                raise ValueError("MD calculations are active, but no configurations read, you may need to wait for calculations to finish." )
        
        if pre_define_configs is not None:
            configs_for_build: list[ase.Atoms] = pre_define_configs
        else:
            configs_for_build: list[ase.Atoms] = self.active_MD_configs
            print(f'Number of active configs = number of QM calculation directories = {len(self.active_MD_configs)}')
        
        for i,config in enumerate(configs_for_build):
            extension = '' if len(configs_for_build)==1 else f'_c_#{i}'
            if QMcode == 'vasp':
                new_calc_dir = vasp.build_vasp_calculation(QM_base_dir,config,f'{label}{extension}',outdir)
            elif QMcode not in __QMcodes__:
                raise ValueError(f'QM code {QMcode} not supported')
            
            self.initialized_QM_dirs.append(new_calc_dir)
            
        print(f'Calculations stored in {outdir}')
        self.QM_base_dir = QM_base_dir
        
        return None
    
    def write_QM_submission_scripts(self, QM_cmd_line: str,
                                   time_per_partition: str,
                                   npartitions: int=1,
                                   save_and_remove: bool=True,
                                   QMcode: str='vasp',
                                   database_file: str=None,
                                   scripts_outdir: str='./QM_scripts',
                                   submit: bool=True,
                                   smart_convergence: bool=False,
                                   expected_motif: np.ndarray=None,
                                   pilot_calculations: bool=True,
                                   mark_as_active: bool=True,
                                   calcs_outdir: str='./QM_calculations',
                                   dependencies: list[str] = []) -> None:
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
        self.hpc_config['dependencies'] = dependencies
        self.hpc_config['time'] = time_per_partition
        header = hpc_utils.fetch_hpc_header(self.hpc_config)
        
        if not self.initialized_QM_dirs:
            self.set_init_QM_dirs(outdir=calcs_outdir)
        #cmd
        if QMcode not in __QMcodes__:
            raise ValueError(f'QM code {QMcode} not supported')
        
        if smart_convergence == True:
            if expected_motif is None:
                raise ValueError('Cannot perform smart convergence without an expected structure (expected_motif).')
            QM_group_indicies, pilot_calculation_configs = group.smart_group_calcs(self.initialized_QM_dirs,ngroups=npartitions,expected_motif=expected_motif,calc_code=QMcode,pilot_calculations=pilot_calculations)
            if pilot_calculations:
                n_main_calcs = len(self.initialized_QM_dirs)
                QM_group_indicies_extended = np.zeros((QM_group_indicies.shape[0],QM_group_indicies.shape[1]+1),dtype=np.int16)
                QM_group_indicies_extended[:,1:] = QM_group_indicies[:,0:]
                for i in range(npartitions):
                    self.build_QM_calculations(self.initialized_QM_dirs[QM_group_indicies[i,0]],QMcode,calcs_outdir,pre_define_configs=[pilot_calculation_configs[i]],label=f'pilot_#{i}')
                    #pilot calculations are the last items in the list. 
                    QM_group_indicies_extended[i,0] = n_main_calcs + i
                QM_group_indicies = QM_group_indicies_extended.copy()

            if QMcode=='vasp':
                for i,dir in enumerate(self.initialized_QM_dirs):
                    if i in list(QM_group_indicies[:,0]):
                        vasp.set_icharg(2,dir)
                    else:
                        vasp.set_icharg(1,dir)
                        
            self.initialized_QM_dirs = [self.initialized_QM_dirs[i] for sublist in QM_group_indicies for i in sublist]
            
        cmd_scipts = write_run_calculation_scripts(self.initialized_QM_dirs,
                                                   QM_cmd_line,
                                                   npartitions=npartitions,
                                                   save_and_remove=save_and_remove,
                                                   python_env=self.hpc_config['python_env'],
                                                   code=QMcode,database_file=database_file,
                                                   smart_convergence=smart_convergence)
        Path(scripts_outdir).mkdir(exist_ok=True)
        for i,cmd in enumerate(cmd_scipts):
            script = hpc_utils.ScriptBuilder(header=header)
            script.add_cmd_line(cmd)
            script.write_script(f'{scripts_outdir}/QM_submission_script_#{i}')
            if submit:
                script.submit_script()
        
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
            self.active_MD_configs.extend(lammps.read_lammps_output(dir,atom_types=self.atom_types))
            
        return self.active_MD_configs
    
    
    def set_active_MD_dirs(self, outdir: str='./MD_calculations') -> None:
        '''
        Set MD calculation directories manually.
        
        Parameters
        ----------
        outdir: str
            path to directory storing MD calculations, expected file structure is {outdir}/{calculation_dir}.
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
        configs = ase.io.read(config_file,':')
        self.active_MD_configs.extend(configs)
        print(f'Number of active configs updated from {old_len} to {len(self.active_MD_configs)}')
        
        return None

    
    def check_initialized_MD_dirs(self):
        '''
        prints initialized MD directories
        '''
        print('MD directories initialized: ')
        print(self.initialized_MD_dirs)
        
    
    def check_active_MD_dirs(self):
        '''
        prints active MD directories
        '''
        print('MD directories currently active: ')
        print(self.active_MD_dirs)
        print('Num active directories: ', len(self.active_MD_dirs))
        
    def check_QM_base_dir(self):
  
        if self.QM_base_dir is None:
            dir = input('Searched for a Quantum mechanical base directory, none found please define a path: ')
            if not Path(dir).is_dir():
                raise FileNotFoundError(f'{dir} not found.')
            else:
                self.QM_base_dir = dir
                return True
        else:
            return True
        

        
        
    def check_active_MD_configs(self):
        '''
        prints active MD configurations.
        '''
        print('MD configurations currently active: ')
        print(self.active_MD_configs)
        print('Num active configurations: ', len(self.active_MD_configs))
        
    def save_active_MD_configs(self, output_file: str):
        '''
        Saves the active MD configurations.
        '''
        
        ase.io.write(output_file,self.active_MD_configs)
        
        
    def set_init_QM_dirs(self, outdir: str='./QM_calculations') -> None:
        '''
        Set QM calculation directories manually.
        '''
        
        path = Path(outdir)
        if not path.exists():
            raise FileNotFoundError(f'Tried to read QM calculations in {outdir} but the directory was not found, may need to set param `calc_dirs`')
        subdirs = [p for p in path.iterdir() if p.is_dir()]
        for calc in subdirs:
            if calc not in self.initialized_QM_dirs:
                self.initialized_QM_dirs.append(str(calc))
                
        return None
          
        
def write_run_calculation_scripts(calc_dirs: list[str],
                                 cmd_line: str,
                                 npartitions: int=1, 
                                 save_and_remove: bool=False,
                                 python_env: str=None,
                                 code: str=None,
                                 database_file: str=None,
                                 smart_convergence: bool=False) -> list[str]:
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
        
        savedata_cmd = f'echo "Saving data from this calculation to {database_file}"\n'
        savedata_cmd += f'{python_env}/bin/python -m mlipts.append_to_database $dir {database_file} {code}'
        remove_cmd = f'echo "Data saved, deleting {code} directory"\n'
        remove_cmd += 'rm -r $dir'
    else:
        savedata_cmd = ''
        remove_cmd = ''
      
    if smart_convergence == True:
        if code == 'vasp':
            smart_convergence_cmd = 'if (( i < num_dirs-1 )); then\necho "Coping CHGCAR to next calculation"\ncp "$dir/CHGCAR" "${directories[i+1]}/"\nfi'
        else:
            raise ValueError(f'smart convergence only supported for [vasp]')
    else:
        smart_convergence_cmd = ''
    
    num_calcs_per_submission = int(len(calc_dirs) / npartitions)

    for i in range(npartitions):
        
        script = ''

        current_dirs=''
        for dir in calc_dirs[int(i*num_calcs_per_submission):int((i+1)*num_calcs_per_submission)]:
            current_dirs+=f'{dir} '
        script+='\n'
        script+=f'directories=({current_dirs})\n'
        script+='num_dirs=${#directories[@]}\n'
        script+='for ((i=0; i<num_dirs; i++)); do\ndir="${directories[i]}"\n'
        script+=f'echo "Running {code} in $dir"\n'
        script+=f'cd $dir\n{cmd_line}\ncd -\n'
        script+=f'{smart_convergence_cmd}\n{savedata_cmd}\n{remove_cmd}\n'
        script+=f'done\n'

        calc_scripts.append(script)
            
    return calc_scripts
            
            
def run_labelling():
    
    print('==============LABELLING CONFIGURATIONS WITH MLIPTS==============')
    
    return None