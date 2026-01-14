'''
@author William Davie

High performance computing utilities.
'''


import yaml
import subprocess
from typing import TypedDict, Optional
from mlipts.hpc_submission import archer2
from mlipts.constants import __hpcs__

__config_keys__ = ['hpc','hpc_account','processor','nodes','ranks','gpus','time','python_env','custom_header','dependencies']

class hpcConfig(TypedDict):
    hpc: str
    hpc_account: str
    processor: str='cpu'
    nodes: int=1
    ranks: int=1
    gpus: int=1
    time: str='01:00:00'
    python_env: str
    custom_header: str = None
    dependencies: list[str]=[]
    

def load_hpc_config(self, hpc_config: str=None, **kwargs) -> hpcConfig:
    '''
    Given a hpc_config.yaml file, use to set attributes of self and return a hpcConfig dict.
    '''

    config_dict = {}
    if hpc_config:
        with open(hpc_config, 'r') as f:
            config_dict = yaml.safe_load(f)
        
    config_dict.update(kwargs)
        
    for key in config_dict:
        if key not in __config_keys__:
            raise ValueError(f'Invalid hpc parameter ({key}) in config file or kwargs, options include: {__config_keys__}')
        setattr(self, key, config_dict[key])
    
    return config_dict
    

def fetch_hpc_header(hpc_config: hpcConfig) -> str:
    '''
    Given some hpc parameters returns the header of a submission script. 
    '''

    hpc = hpc_config['hpc']
    if hpc == 'archer2':
        if hpc_config['processor'] == 'cpu':
            header = archer2.archer2_submission_template(nodes=hpc_config['nodes'],ranks=hpc_config['ranks'],time=hpc_config['time'],account=hpc_config['hpc_account'],dependencies=hpc_config['dependencies'])
        elif hpc_config['processor'] == 'gpu':
            header = archer2.archer2_gpu_submission_template(account=hpc_config['hpc_account'],time=hpc_config['time'],gpus=hpc_config['gpus'],dependencies=hpc_config['dependencies'])
        else:
            raise ValueError(f'Cannot accept processor named {hpc_config["processor"]}, choose "cpu" or "gpu".')
    
    elif hpc == 'custom':
        if hpc_config['custom_header'] == None:
            raise ValueError('custom hpc header but no header_str argument provided.')
        else:
            header = hpc_config['custom_header']
    
    elif hpc not in __hpcs__:
        raise ValueError(f'hpc {hpc} not supported.')
    
    return header

'''
A class for building hpc submission scripts.
'''

class ScriptBuilder():
    
    def __init__(self, header) -> None:
        self.script_lines = []
        self.output_path: str = None
        self.header = header
        return None
    
    def add_cmd_line(self, cmd_line: str):
        '''
        adds line to script object
        '''
        self.script_lines.append(cmd_line)
        return None
    
    def write_script(self, output_path: str, message: bool=True):
        if self.header == '':
            print("<!> Warning : writing a script with no header.")
        if not self.script_lines:
            print('<!> Warning : writing an empty script')

        with open(output_path,'w') as f:
            f.write(self.header)
            for line in self.script_lines:
                f.write(line + '\n')
                
        self.output_path = output_path
        if message:  
            print(f'Job submission script saved to {output_path}.')
        
        return None
    
    def submit_script(self) -> int:
        '''
        submits the script with sbatch
        '''
        result = subprocess.run(f'sbatch {self.output_path}',shell=True,capture_output=True,text=True)
        print(result.stdout)
        try: return int(result.stdout.split()[-1])
        except:
            print('<!> Warning slurm job id not found')
            return None
    
    def add_loop():
        '''
        Add a loop like:
            script.add_cmd_line(f'configs=({configs_to_run})')
            script.add_cmd_line('num_configs=${#configs[@]}')
            script.add_cmd_line('for ((i=0; i<num_configs; i++)); do\nconfig="${configs[i]}"')
            script.add_cmd_line('echo "Training model with configuration $config"')
            script.add_cmd_line(cmd_line + ' $config')
            script.add_cmd_line(f'done')
        '''
        return None