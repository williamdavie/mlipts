'''

@author William Davie

With a calculation done, we now want to append positions, energies and forces to the training database.
Formated for MACE, details on data format:  https://github.com/ilyes319/mace-tutorials
'''

import numpy as np
import sys,os
from mlipts.constants import __MDcodes__,__QMcodes__
from mlipts.codes.vasp import append_vasp_calc_to_database

def append_to_database(database_file: str, calc_dir: str, code: str='vasp', forces: bool=True, pbc: str='T T T'):
    '''
    Given a calculation in directory calc_dir, read positions and forces and 
    '''
    
    if code == 'vasp':
        append_vasp_calc_to_database(database_file,calc_dir,pbc=pbc)
    else:
        raise ValueError(f'code {code} not supported.')
    
    return None

# called by hpc
if __name__ == '__main__': 
    
    calc_dir = sys.argv[1]
    database_file = sys.argv[2]
    code = sys.argv[3]
    
    append_to_database(database_file,calc_dir,code)