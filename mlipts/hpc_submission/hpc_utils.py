
from mlipts.hpc_submission import archer2
from mlipts.constants import __hpcs__

def fetch_hpc_header(hpc: str, hpc_account: str,
                     processor: str='cpu',
                     nodes: int=1,
                     ranks: int=1,
                     gpus: int=1,
                     time: str='01:00:00',
                     header_str: str=None) -> str:
    '''
    Given some hpc parameters returns the header of a submission script. 
    '''

    if hpc == 'archer2':
        if processor == 'cpu':
            header = archer2.archer2_submission_template(nodes=nodes,ranks=ranks,time=time,account=hpc_account)
        elif processor == 'gpu':
            header = archer2.archer2_gpu_submission_template(account=hpc_account,time=time,gpus=gpus)
        else:
            raise ValueError(f'Cannot accept processor named {processor}, choose "cpu" or "gpu".')
    
    elif hpc == 'custom':
        if header_str == None:
            raise ValueError('custom hpc header but no header_str argument provided.')
        else:
            header = header_str
    
    elif hpc not in __hpcs__:
        raise ValueError(f'hpc {hpc} not supported.')
    
    return header