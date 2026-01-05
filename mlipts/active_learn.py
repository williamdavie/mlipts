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
from mlipts.hpc_submission import hpc_utils
from mlipts.constants import __architectures__

class ActiveLearn():
    
    def __init__(self):
        
        self.active_model_configs = []
        
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
            
        Path(outdir).mkdir(exist_ok=True)
        for i in range(n_models):
            base_config['name'] = f'model_#{i}'
            base_config['seed'] = np.random.randint(0,1000)
            new_config_file = outdir+f'/config_#{i}'
            
            with open(new_config_file,'w') as f:
                yaml.safe_dump(base_config,f)
                
            self.active_model_configs.append(new_config_file)
            
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
                       outdir: str='training_scripts'):
        '''
        Generates submission scripts to train the set of models. 
        '''
        
        Path(outdir).mkdir(exist_ok=True)
        
        if not self.active_model_configs:
            self.fetch_model_configs()
        
        # header
        header = hpc_utils.fetch_hpc_header(hpc=hpc,hpc_account=hpc_account,processor=processor,nodes=nodes,ranks=ranks,gpus=gpus,time=time,header_str=custom_header_str)
        
        if architecture == 'mace':
            header += '\n'
            header += f'source {python_env}/bin/activate\n' # mace is written in python. 
            cmd_line = f'{python_env}/bin/mace_run_train --config'
        else:
            raise ValueError(f'Recieved unkown architechture {architecture}, currently availible architechtures are {__architectures__}')
        
        num_calcs_per_submission = int(len(self.active_model_configs) / npartitions)
        
        for i in range(npartitions):
            script = ''
            configs_to_run =''
            for config in self.active_model_configs[int(i*num_calcs_per_submission):int((i+1)*num_calcs_per_submission)]:
                configs_to_run+=f'{config} '
                
            script+=f'configs=({configs_to_run})\n'
            script+='num_configs=${#configs[@]}\n'
            script+='for ((i=0; i<num_configs; i++)); do\nconfig="${configs[i]}"\n'
            script+='echo "Training model with configuration $config"\n'
            script+=cmd_line + ' $config\n'
                
            with open(f'{outdir}/train_script_#{i}','w') as f:
                f.write(header)
                f.write(script)
                
            print(f'Submission script saved to {outdir}/train_script_#{i}.')
            
            if submit == True:
                subprocess.run(f'sbatch {outdir}/train_script_#{i}',shell=True)
            else:
                print(f'You opted not to batch the submission script, it can batched via the working directory with cmd line: \n sbatch {outdir}/train_script_#{i}')

        return None
    
    
    
    # generates a set of lammps calculations to use the newly trained models
    def build_commitee_test_suite():
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