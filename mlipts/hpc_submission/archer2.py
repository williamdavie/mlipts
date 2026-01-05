import re

def archer2_submission_template(nodes: int, ranks: int, time: str, account: str):
    
    if not re.match(r'^\d{2}:\d{2}:\d{2}$', time):
        raise ValueError('Time must have format XX:XX:XX')
    
    hours, minutes, seconds = time.split(":")
    
    if int(minutes) <= 20 and int(hours) == 0:
        
        qos='short'

    else:
        
        qos = 'standard'
    
    
    return f'''#!/bin/bash

#SBATCH --job-name=job_MLIPTS
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={ranks}
#SBATCH --cpus-per-task=1
#SBATCH --time={time}

#SBATCH --account={account}
#SBATCH --partition=standard
#SBATCH --qos={qos}
    '''


def archer2_gpu_submission_template(account: str, time: str, nodes: int=1, gpus: int=1, load_default_modules: bool=True):
    
    if load_default_modules:
        module_load_str = '''module load PrgEnv-amd
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-milan'''
    
    else:
        module_load_str = ''
        
    return f'''#!/bin/bash

#SBATCH --job-name=MLIPTS_gpu
#SBATCH --account={account}
#SBATCH --partition=gpu
#SBATCH --qos=gpu-shd
#SBATCH --nodes={nodes}
#SBATCH --gpus={gpus}
#SBATCH --time={time}

{module_load_str}
    '''
    
    