'''
@author William Davie

filter functions 
future improvements include additional filter methods and optimization.
'''

from ase import Atoms
import numpy as np
from mlipts.similarity.emd import EMD
from mlipts.similarity.pdd import atoms_configs_PDDs
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt

def filter_by_emd(configs: list[Atoms], tol: float, k: int, show_dendrograms: bool=False) -> tuple:
    '''
    filters a set of configurations by earth movers distance according to a tolerance 
    '''
    
    PDDs = atoms_configs_PDDs(configs,k)
    emds = []
    emds_filtered = []
    
    remove_ = set()
    for i in range(len(PDDs)):
        print(f"\rProgress: {int(100*i/len(PDDs))}%", end="")
        if i in remove_:
            continue
        for j in range(i+1, len(PDDs)):
            emd = EMD(PDDs[i], PDDs[j])
            if j in remove_:
                continue
            if emd <= tol:
                remove_.add(j)
            else:
                emds_filtered.append(emd)          
    print('\n')

    inds = [i for i in range(len(PDDs)) if i not in remove_]
    
    return [configs[i] for i in inds], inds


