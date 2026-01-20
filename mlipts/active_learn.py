'''
@author William Davie

Main class for the active learning process. Currently MACE specific. 

Dev note:
There are a number of parallels between DataCollection and ActiveLearn. May be possible to define a parent class
'''

import yaml
import numpy as np
from pathlib import Path
import subprocess
import ase
import ase.io
from mlipts.hpc_submission import hpc_utils
from mlipts.data_collection import DataCollection
from mlipts.constants import __architectures__

class ActiveLearn():
    
    def __init__(self, hpc_config: str=None, architecture: str='mace', **kwargs):
        
        self.active_model_config_files: list[str] = []
        self.models_dir: str = None
        self.hpc_config = hpc_utils.load_hpc_config(self,hpc_config, **kwargs)
        self.architecture = architecture
        
        if self.architecture == 'mace':
            if not self.hpc_config['python_env']:
                raise ValueError('To use mace for active learning you must define a python env with mace installed.')
        elif self.architecture not in __architectures__:
            raise ValueError(f'Active learning for architechture {architecture} not supported, supported architectures: {__architectures__}') 
        
        self.expected_final_models = []
        
        return None
    
    def define_commitee(self, base_config_file: str, n_models: int, outdir: str='model_configs') -> None:
        '''
        Given a model configuration file and the number of desired models, generates a commitee config files, randomizing the seed.
        
        Parameters
        ----------
        base_config_file: str
            filename of a config.yml file, the settings of this file are copied to generate a model comittee.
        n_models: int
            number of models you want to generate.
            
        Returns
        -------
        None: None
            generates comittee of config.yml files. 
        '''
        
        with open(base_config_file,'r') as f:
            base_config = yaml.safe_load(f)
            
        self.models_dir = base_config['model_dir']
            
        Path(outdir).mkdir(exist_ok=True)
        for i in range(n_models):
            base_config['name'] = f'model_{i}'
            base_config['seed'] = np.random.randint(0,1000)
            new_config_file = outdir+f'/config_#{i}'
            
            with open(new_config_file,'w') as f:
                yaml.safe_dump(base_config,f)
                
            self.active_model_config_files.append(new_config_file)
            
        self.expected_final_models = [f'model_{i}_stagetwo.model' for i in range(n_models)]
            
        return None
    
    # generates submission scipts to train the model
    def train_commitee(self,
                       time_per_partition: str,
                       npartitions: int=1,
                       submit: bool=True,
                       outdir: str='scripts',
                       dependencies: list[str]=[]):
        '''
        Generates submission scripts to train the set of models. 
        '''
        
        Path(outdir).mkdir(exist_ok=True)
        
        if not self.active_model_config_files:
            self.fetch_model_configs()
        
        # header
        self.hpc_config['time'] = time_per_partition
        self.hpc_config['dependencies'] = dependencies
        header = hpc_utils.fetch_hpc_header(self.hpc_config)
        
        if self.architecture == 'mace':
            header += '\n'
            header += f'source {self.hpc_config["python_env"]}/bin/activate\n' # mace is written in python. 
            cmd_line = f'{self.hpc_config["python_env"]}/bin/mace_run_train --config'
        
        
        num_calcs_per_submission = int(len(self.active_model_config_files) / npartitions)
        
        for i in range(npartitions):
            script = hpc_utils.ScriptBuilder(header=header)
            configs_to_run =''
            for config in self.active_model_config_files[int(i*num_calcs_per_submission):int((i+1)*num_calcs_per_submission)]:
                configs_to_run+=f'{config} '
                
            script.add_cmd_line(f'configs=({configs_to_run})')
            script.add_cmd_line('num_configs=${#configs[@]}')
            script.add_cmd_line('for ((i=0; i<num_configs; i++)); do\nconfig="${configs[i]}"')
            script.add_cmd_line('echo "Training model with configuration $config"')
            script.add_cmd_line(cmd_line + ' $config')
            script.add_cmd_line(f'done')
                
            script.write_script(f'{outdir}/train_script_#{i}')
        
            if submit == True:
                script.submit_script()
                

        return None
    
    def evaluate_committee(self,
                           time_per_partition: str,
                           npartitions: int=1,
                           atomic_config_file: str='MD_samples.xyz', # consistency with data_collection class.
                           submit: bool=True,
                           models_dir: str=None, script_outdir='scripts',
                           dependencies: list[str]=[],
                           create_lammps_model: bool=True):
        '''
        Given a config file, evaluates the energy and forces on each configuration using each model in the committee. 
        '''
        
        Path(f'{script_outdir}').mkdir(exist_ok=True)
        Path('evaluated_configs').mkdir(exist_ok=True)
            
        self.hpc_config['time'] = time_per_partition
        self.hpc_config['dependencies'] = dependencies
        header = hpc_utils.fetch_hpc_header(self.hpc_config)
        
        if (self.models_dir is None and models_dir is None) or not Path(models_dir).exists():
            raise ValueError('You asked to evaluate a committee models but have not specified where the models are located')
        elif models_dir is not None:
            self.models_dir = models_dir 
    

        all_model_files = [str(p) for p in Path(self.models_dir).iterdir()]
        
        if self.architecture == 'mace':
            header += '\n'
            header += f'source {self.hpc_config["python_env"]}/bin/activate\n'
            
            
            if self.expected_final_models:
                models = [f'{self.models_dir}/{self.expected_final_models[i]}' for i in range(len(self.expected_final_models))]
            else:
                models = [i for i in all_model_files if ('stagetwo.model' in i and 'run' not in i)]
  
            eval_model_cmd = f'{self.hpc_config["python_env"]}/bin/mace_eval_configs --configs {atomic_config_file} --model $model --output $model_output'
        
        if create_lammps_model:
            lmps_model_cmd = 'echo "Creating a lammps model"\npython -m mace.cli.create_lammps_model $model'
        else:
            lmps_model_cmd = ''
        
        num_calcs_per_submission = int(len(models) / npartitions)
        
        
        for i in range(npartitions):
            script = hpc_utils.ScriptBuilder(header=header)
            models_to_evaluate = ''
            for model in models[int(i*num_calcs_per_submission):int((i+1)*num_calcs_per_submission)]:
                models_to_evaluate+=f'{model} '
                
            script.add_cmd_line(f'models=({models_to_evaluate})')
            script.add_cmd_line('num_models=${#models[@]}')
            script.add_cmd_line('for ((i=0; i<num_models; i++)); do')
            script.add_cmd_line('model="${models[i]}"\nmodel_output="evaluated_configs/configs_$i.xyz"')
            script.add_cmd_line('echo "Evaluating $model"')
            script.add_cmd_line(eval_model_cmd)
            script.add_cmd_line(lmps_model_cmd)
            script.add_cmd_line('done')
            
            script.write_script(f'{script_outdir}/evaluate_script_#{i}')
            
            if submit == True:
                script.submit_script()

    def uncertainty_script(self,tol,time: str='00:20:00',dependacy=None):
        '''
        Writes a short script to quantify uncertainty on an hpc.
        '''
        self.hpc_config['time'] = time
        header = hpc_utils.fetch_hpc_header(self.hpc_config)
        
    
    # generates the submission scripts for above. saves the data from each model to xyz
    def test_suite_submission():
        
        return None
    
    
    # filters the configs from each model and defines a final training set. (Then data class can be used)
    def run_train():
        return None
    

    def fetch_model_configs():
        return None
    

        
    

class UncertaintyQuantification():
    
    def __init__(self, original_sample_file: str, evaluated_config_files: list[str], architecture='mace'):
        
        self.original_sample_file = original_sample_file
        self.evaluated_config_files = evaluated_config_files
        self.committee_size = len(evaluated_config_files)
        self.architecture = architecture
        self.uncertainties = None
        
        if architecture == 'mace':
            self.energy_tag = 'MACE_energy'
            self.forces_tag = 'MACE_forces'
        else:
            raise ValueError(f'Architechture {architecture} not supported')
        
        self.original_configs =  ase.io.read(original_sample_file,':')
        self.n_configs = len(self.original_configs)
        self.all_configs = np.empty((self.n_configs,self.committee_size),dtype=ase.Atoms)
        
        for i in range(self.committee_size):
            current_configs = ase.io.read(evaluated_config_files[i],':')
            
            if len(current_configs) != self.all_configs.shape[0]: 
                raise ValueError('Files in the committee have a different number of evaluated configurations')
            
            self.all_configs[:,i] = current_configs
            
    def filter_configs(self, tol: float, max_configs: int=None, min_configs: int=None,method='dubois',selection='largest') -> None:
        
        if self.uncertainties is None:
            if method == 'dubois':
                self.uncertainties = self.dubois_uncertainty()
            else:
                raise ValueError(f'Method: {method}, unknown')
        
        self.indices = np.where(self.uncertainties > tol)[0]
        if max_configs is not None:
            if len(self.indices) > max_configs:
                self.indices = np.argsort(self.uncertainties)[::-1][0:max_configs]
            
        if min_configs is not None:
            if len(self.indices) < min_configs:
                self.indices = np.argsort(self.uncertainties)[::-1][0:min_configs]
        
        new_configs = [self.original_configs[i] for i in self.indices]
        
        ase.io.write(f'active_learning_result.xyz',new_configs)
        
        return None
        
            
    def dubois_uncertainty(self) -> tuple:
    
        energy_deviations = np.zeros(self.n_configs)
        force_deviations = np.zeros(self.n_configs)
        
        for i in range(self.n_configs):
            current_configs = self.all_configs[i]
            energy_deviations[i] = (self.energy_standard_deviation(current_configs))
            force_deviations[i] = (self.ave_force_standard_deviation(current_configs))
          
        all_energy_sd = np.std(energy_deviations)
        all_force_sd = np.std(force_deviations)
        
        dubois_uncertainties = energy_deviations / all_energy_sd + force_deviations / all_force_sd
        
        self.uncertainties = dubois_uncertainties
        
        return dubois_uncertainties


    def energy_standard_deviation(self,configs: list[ase.Atoms]) -> float:
        '''
        Given a list of equal configurations with different model energies, return a standard deviation.
        '''
        energies = np.zeros(self.committee_size)
        for i in range(self.committee_size): energies[i] = configs[i].info[self.energy_tag]
        
        return np.std(energies)
    

    def ave_force_standard_deviation(self,configs: list[ase.Atoms]) -> float:
        '''
        Given a list of equal configurations with different forces, return a standard deviation.
        '''
        n_atoms = len(configs[0])
        forces = np.zeros((self.committee_size,n_atoms,3))
        
        for i in range(self.committee_size): forces[i] = configs[i].arrays[self.forces_tag]
        
        sds = np.zeros((n_atoms,3)) # list of standard deviations
        for i in range(n_atoms):
            for j in range(3):
                sd = np.std(forces[:,i,j])
                sds[i,j] = sd
        
        return np.average(sds)

    


def run_active_learn(**kwargs) -> list[ase.Atoms]:
    '''
    Given a set of base_config, training, test data, runs the full active learning workflow from start to finish. The result is a set of configurations to be labelled. 
    
    Steps
    -----
    1. A commitee of configs are generated.
    2. A commitee of models are trained.
    3. A set of new configurations are sampled.
    4. New configurations are evaluated by the committee.
    5. Uncertainty is quantified and the new configurations are given as output.
    '''       

    print('==============ACTIVE LEARNING WITH MLIPTS==============')
    print('\n')
    print('Architechture being used: MACE')
    print('\n')
    hpc_config = input('Input name of HPC configuration file: ')
    activelearning = ActiveLearn(hpc_config,'mace',**kwargs)
    print('------------------Committee setup---------------------')
    committee_size = input('Input committee size: ')
    base_config = input('Input name of configuration file: ')
    activelearning.define_commitee(base_config,committee_size)
    print(f'Copying {base_config} {committee_size} times with a random seed.')
    print('\n')
    print('Setting up a script to begin model training.')
    training_job_id = 0 
    print(f'Training at job id: ')
    print('------------Sampling new configurations---------------------')
    defined = input('Have you already sampled a new configurations ("y"/"n") ')
    if defined.capitalize() == 'Y':
        new_configs_file = input('Input location of new configurations file: ')
    else: 
        sampling_job_id = 0
    print('------------Evaluate models and quantify uncertainty--------------------')
    output_name = input('Input')
    print(f'Setting up scripts to evaluate models with dependacies: {training_job_id, sampling_job_id}')
    
    return None
    