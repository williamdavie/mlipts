'''

@author William Davie

'''

from ase import Atoms
import numpy as np
from mlipts.similarity.emd import EMD
from mlipts.similarity.pdd import atoms_configs_PDDs


def filter_by_emd(configs: list[Atoms], tol: float, k: int) -> tuple:
    '''
    filters a set of configurations by earth movers distance according to a tolerance 
    '''
    
    PDDs = atoms_configs_PDDs(configs,k)
    emds = []
    
    remove_ = set()
    for i in range(len(PDDs)):
        if i in remove_:
            continue
        for j in range(i+1, len(PDDs)):
            if j in remove_:
                continue
            emd = EMD(PDDs[i], PDDs[j])
            if emd <= tol:
                remove_.add(j)
    
    inds = [i for i in range(len(PDDs)) if i not in remove_]
    
    return [configs[i] for i in inds], inds


