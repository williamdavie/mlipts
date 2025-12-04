
'''
@author William Davie

Main class for collecting a data set for model training/fine-tuning. 
'''

from mlipts.codes.lammps import build_lammps_calculations, write_lammps_bash_command, read_lammps_output
from mlipts.hpc_submission.archer2 import archer2_submission_template
from mlipts.similarity.filter import filter_by_emd
from ase import Atoms
from ase.io import read, write

__hpcs__ = ['Archer2']
__MDcodes__ = ['lammps']
__QMcodes__ = ['VASP']
__diffmethods__ = ['emd']


class DataCollection():
    
    
    def __init__(self, atom_types: list[str]):
        '''
        Initialize a DataCollection class. Used to develop a database for training a Machine Learned Interatomic Potential.
        '''

        self.atom_types = atom_types
        
        self.complete_MD = []
        self.active_MD = []
        self.initialized_MD = []
        self.MD_submission_count = 0
        
        self.active_MD_configs = []
        
    def build_MD_calculations(self, MD_base_dir: str, MDcode: str='lammps', outdir: str='.') -> None:
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
            new_dirs = build_lammps_calculations(MD_base_dir,None,outdir=outdir)
            self.initialized_MD.extend(new_dirs)
        
        elif MDcode not in __MDcodes__:
            raise ValueError(f'{MDcode} not supported.')
        
        return None
    
    def build_MD_submission_script(self, MD_cmd_line: str,
                                   nodes: int, ranks: int,
                                   time: str,
                                   MDcode: str='lammps',
                                   hpc: str='Archer2',
                                   mark_as_active: bool=True):
        '''
        Build submission script for Molecular Dynamics simulations, built for all directories marked 'incomplete'. 
        
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
        
        Returns
        -------
        None : None
            generates MD_submission_script_#i in the working directory.
        
        '''
        # header
        if hpc == 'Archer2':
            header = archer2_submission_template(nodes,ranks,time)
            
        elif hpc not in __hpcs__:
            raise ValueError(f'hpc {hpc} not supported.')
        
        # cmd
        if MDcode == 'lammps':
            cmd = write_lammps_bash_command(self.initialized_MD,MD_cmd_line)
    
        elif MDcode not in __MDcodes__:
            raise ValueError(f'MD code {MDcode} not supported.')
        

        with open (f'MD_submission_script_#{self.MD_submission_count}','w') as f:
            
            f.write(header)
            f.write('\n')
            f.write(cmd)
            
            
        print(f'MD submission script saved to: MD_submission_script_#{self.MD_submission_count}')
        self.MD_submission_count += 1
        
        if mark_as_active:
            self.active_MD.extend(self.initialized_MD)
            self.initialized_MD = []
      
        return None
    

    def filter_active_MD(self,tol: float, method: str='emd', auto: bool=False, k: int=20, show_dendrograms: bool=False):
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
        
        self.fetch_active_MD_configs()
        
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
    
    
    def fetch_active_MD_configs(self):
        '''
        Collect the atomic configs from the active MD directories
        '''
        self.active_MD_configs = []

        for dir in self.active_MD:
            self.active_MD_configs.extend(read_lammps_output(dir,atom_types=self.atom_types))
        return self.active_MD_configs


    def mark_MD_complete():
        return None
        
    
    def build_QM_calculations(all: bool=False):
       
        return None
    
    def build_QM_submission_scripts():
        
        return None
    
    def __check_initialized__(self):
        '''
        prints initialized MD directories
        '''
        print('MD directories initialized: ')
        print(self.initialized_MD)
        
    
    def __check_active__(self):
        '''
        prints active MD directories
        '''
        print('MD directories currently active: ')
        print(self.active_MD)
        
    def __save_active__(self, outname: str):
        
        write(f'{outname}.xyz',self.active_MD_configs)
        
        