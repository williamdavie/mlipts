import re

def archer2_submission_template(nodes: int, ranks: int, time: str, account: str, dependencies: list[str]=[]):
    
    if not re.match(r'^\d{2}:\d{2}:\d{2}$', time):
        raise ValueError('Time must have format XX:XX:XX')
    
    hours, minutes, seconds = time.split(":")
    
    if int(minutes) <= 20 and int(hours) == 0:
        
        qos='short'

    else:
        
        qos = 'standard'
        
    if dependencies:
        dependency_str = '#SBATCH --dependency=afterok'
        for id in dependencies:
            dependency_str += f':{id}'
    else:
        dependency_str = '\n'
    
    
    return f'''#!/bin/bash

#SBATCH --job-name=job_MLIPTS
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={ranks}
#SBATCH --cpus-per-task=1
#SBATCH --time={time}
#SBATCH --account={account}
#SBATCH --partition=standard
#SBATCH --qos={qos}
{dependency_str}
'''


def archer2_gpu_submission_template(account: str, time: str, nodes: int=1, gpus: int=1, dependencies: list[str]=[],load_default_modules: bool=True):
    
    if load_default_modules:
        module_load_str = '''module load PrgEnv-amd
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-milan'''
    
    else:
        module_load_str = ''
        
    if dependencies:
        dependency_str = '#SBATCH --dependency=afterok'
        for id in dependencies:
            dependency_str += f':{id}'
    else:
        dependency_str = ''
        
    return f'''#!/bin/bash

#SBATCH --job-name=MLIPTS_gpu
#SBATCH --account={account}
#SBATCH --partition=gpu
#SBATCH --qos=gpu-shd
#SBATCH --nodes={nodes}
#SBATCH --gpus={gpus}
#SBATCH --time={time}
{dependency_str}
{module_load_str}
    '''
    
    