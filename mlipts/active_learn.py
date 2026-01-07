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
from mlipts.constants import __architectures__

class ActiveLearn():
    
    def __init__(self):
        
        self.active_model_config_files: list[str] = []
        self.models_dir: str = None
        
        return None
    
    
    def define_commitee(self, base_config_file: str, n_models: int, outdir: str='model_configs',architecture: str='MACE') -> None:
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
            base_config['name'] = f'model_#{i}'
            base_config['seed'] = np.random.randint(0,1000)
            new_config_file = outdir+f'/config_#{i}'
            
            with open(new_config_file,'w') as f:
                yaml.safe_dump(base_config,f)
                
            self.active_model_config_files.append(new_config_file)
            
        return None
    
    # generates submission scipts to train the model
    def train_commitee(self, hpc: str,
                       hpc_account: str,
                       time: str,
                       processor: str='gpu',
                       nodes: int=1,
                       ranks: int=1, # if cpu selected
                       gpus: int=1, # if gpu selected
                       npartitions: int=1,
                       python_env: str=None,
                       architecture: str='mace',
                       submit: bool=True,
                       custom_header_str: str=None,
                       outdir: str='scripts'):
        '''
        Generates submission scripts to train the set of models. 
        '''
        
        Path(outdir).mkdir(exist_ok=True)
        
        if not self.active_model_config_files:
            self.fetch_model_configs()
        
        # header
        header = hpc_utils.fetch_hpc_header(hpc=hpc,hpc_account=hpc_account,processor=processor,nodes=nodes,ranks=ranks,gpus=gpus,time=time,header_str=custom_header_str)
        
        if architecture == 'mace':
            header += '\n'
            header += f'source {python_env}/bin/activate\n' # mace is written in python. 
            cmd_line = f'{python_env}/bin/mace_run_train --config'
        else:
            raise ValueError(f'Recieved unkown architechture {architecture}, currently availible architechtures are {__architectures__}')
        
        num_calcs_per_submission = int(len(self.active_model_config_files) / npartitions)
        
        for i in range(npartitions):
            script = ''
            configs_to_run =''
            for config in self.active_model_config_files[int(i*num_calcs_per_submission):int((i+1)*num_calcs_per_submission)]:
                configs_to_run+=f'{config} '
                
            script+=f'configs=({configs_to_run})\n'
            script+='num_configs=${#configs[@]}\n'
            script+='for ((i=0; i<num_configs; i++)); do\nconfig="${configs[i]}"\n'
            script+='echo "Training model with configuration $config"\n'
            script+=cmd_line + ' $config\n'
            script+=f'done\n'
                
            with open(f'{outdir}/train_script_#{i}','w') as f:
                f.write(header)
                f.write(script)
                
            print(f'Submission script saved to {outdir}/train_script_#{i}.')
            
            if submit == True:
                subprocess.run(f'sbatch {outdir}/train_script_#{i}',shell=True)
            else:
                print(f'You opted not to batch the submission script, it can batched via the working directory with cmd line: \n sbatch {outdir}/train_script_#{i}')

        return None
    
    def evaluate_committee(self,config_file: str,
                           hpc_account: str,
                           time: str,
                           processor: str='gpu',
                           nodes: int=1,
                           ranks: int=1, # if cpu selected
                           gpus: int=1, # if gpu selected
                           npartitions: int=1,
                           python_env: str=None,
                           architecture: str='mace',
                           hpc: str='archer2',
                           submit: bool=True,
                           custom_header_str: str=None,
                           models_dir: str=None, evaluted_outdir: str='evaluated_samples', script_outdir='scripts'):
        '''
        Given a config file, evaluates the energy and forces on each configuration using each model in the committee. 
        '''
        
        header = hpc_utils.fetch_hpc_header(hpc=hpc,hpc_account=hpc_account,processor=processor,nodes=nodes,ranks=ranks,gpus=gpus,time=time,header_str=custom_header_str)
        
        if (self.models_dir is None and models_dir is None) or not Path(models_dir).exists():
            raise ValueError('You asked to evaluate a committee models but have not specified where the models are located')
        elif models_dir is not None:
            self.models_dir = models_dir 
    
        all_model_files = [str(p) for p in Path(self.models_dir).iterdir()]
        
        if architecture == 'mace':
            if python_env is None:
                raise ValueError('Mace evaluation requires specification of a python enviroment with mace and mlipts installed.')
            else:
                header += '\n'
                header += f'source {python_env}/bin/activate\n'
            
            models = [i for i in all_model_files if ('stagetwo.model' in i and 'run' not in i)]
            eval_model_cmd = f'{python_env}/bin/mace_eval_configs --configs {config_file} --model $model --output $model_output\n'
        
        else:
            raise ValueError(f'Architechture {architecture} not supported')
        
        num_calcs_per_submission = int(len(models) / npartitions)
        
        # header
       
        # script
        
        for i in range(npartitions):
            script = ''
            models_to_evaluate = ''
            for model in models[int(i*num_calcs_per_submission):int((i+1)*num_calcs_per_submission)]:
                models_to_evaluate+=f'{model} '
                
            script+=f'models=({models_to_evaluate})\n'
            script+='num_models=${#models[@]}\n'
            script+='for ((i=0; i<num_models; i++)); do\n'
            script+='model="${models[i]}"\nmodel_output="evaluated_configs_$i.xyz"\n'
            script+='echo "Evaluating $model"\n'
            script+=eval_model_cmd
            script+=f'done\n'
            
            output_str = f'{script_outdir}/evaluate_script_#{i}'
            
            with open(output_str,'w') as f:
                f.write(header)
                f.write(script)
                
                print(f'Submission script saved to {output_str}.')
            
            if submit == True:
                subprocess.run(f'sbatch {output_str}',shell=True)
            else:
                print(f'You opted not to batch the submission script, it can batched via the working directory with cmd line: \n sbatch {output_str}')

        
                
            
def run_active_learn(hpc, hpc_account) -> list[ase.Atoms]:
    '''
    Given a set of base_config, training, test data, runs the full active learning workflow from start to finish. The result is a set of configurations to be labelled. 
    
    Steps
    -----
    1. A commitee of configs are generated.
    2. A commitee of models are trained.
    3. A set of new configurations are sampled and evaluated by the committee.
    4. Uncertainty is quantified and the new configurations are given as output.
    '''       

    print('==============ACTIVE LEARNING WITH MLIPTS==============')
    print('\n')
    print('------------------Committee setup---------------------')
    committee_size = input('Input committee size: ')
    base_config = input('Input a configuration file: ')
    print('Copying {base_config} {commitee_size} times with a random seed.')
    print('Setting up a script to begin model training.')
    print('------------Sampling new configurations---------------------')
    
    
    return None
        
        

        
            
            
            
        
            
        
        
    
    # generates the submission scripts for above. saves the data from each model to xyz
    def test_suite_submission():
        
        return None
    
    
    # filters the configs from each model and defines a final training set. (Then data class can be used)
    def quantify_uncertainty():
        return None
    
    
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
            
    def filter_configs(self, tol: float, method='dubois') -> None:
        
        if method == 'dubois':
            uncertainties = self.dubois_uncertainty()
        else:
            raise ValueError(f'Method: {method}, unknown')
        
        indices = np.where(uncertainties > tol)[0]
        
        new_configs = [self.original_configs[i] for i in indices]
        
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
            
        return dubois_uncertainties


    def energy_standard_deviation(self,configs: list[ase.Atoms]) -> float:
        '''
        Given a list of equal configurations with different model energies, return a standard deviation.
        '''
        energies = np.zeros(self.committee_size)
        for i in range(self.committee_size): energies[i] = configs[i].arrays[self.energy_tag]
        
        return np.std(energies)
    

    def ave_force_standard_deviation(self,configs: list[ase.Atoms]) -> float:
        '''
        Given a list of equal configurations with different forces, return a standard deviation.
        '''
        n_atoms = configs[0].get_number_of_atoms()
        forces = np.zeros((self.committee_size,n_atoms,3))
        
        for i in range(self.committee_size): forces[i] = configs[i].arrays[self.forces_tag]
        
        sds = np.zeros((n_atoms,3)) # list of standard deviations
        for i in n_atoms:
            for j in range(3):
                sd = np.std(forces[:,i,j])
                sds[i,j] = sd
        
        return np.average(sds)

        
    