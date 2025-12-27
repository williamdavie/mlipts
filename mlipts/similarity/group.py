'''
group atomic positions according to similarity criteria.
'''

import numpy as np
from mlipts.codes.vasp import fetch_configs_vasp
from ase.atoms import Atoms
from itertools import product
from ase.io import read,write
from mlipts.similarity.pdd import PDD
from mlipts.similarity.emd import EMD

def smart_group_calcs(calc_dirs: list[str], ngroups: int, expected_motif: np.ndarray, calc_code: str='vasp', group_by: str='emd'):
    '''
    Given a set of calculation paths (calc_dirs), group to maximise calculation convergence when batched to a supercomputer.
    
    Parameters
    ----------
    calc_dirs: list[str]
        list of paths to each calculation.
    ngroups: int
        number of groups to sort calculations into
    expected_motif: int
        a (N,3) array containin the expected atomic motif, i.e. what a relaxed structure would look like.
    calc_code: str
        code used to perform QM calculations. Default is vasp
    group_by:
        method to measure similarity between configurations and perform grouping. 
        
    Returns
    -------
    calc_dirs_grouped: list[list[str]]
        list of paths, each sub list is a group. 
    '''
    
    if calc_code == 'vasp':
        configs = fetch_configs_vasp(calc_dirs)
    else:
        raise ValueError(f'Calculation code {calc_code} not supported in smart grouping. ')
    if group_by == 'emd':
        k = int(input('Input number of neighbours (k) used for earth movers distance: '))
        group_indicies = smart_group_by_emd(configs,ngroups,expected_motif,k)
    else:
        raise ValueError(f'similarity assessment statergy (group_by) {group_by} not regonised')

    return [calc_dirs[i] for sublist in group_indicies for i in sublist]



def smart_group_by_emd(configs: list[Atoms], ngroups: int, expected_motif: np.ndarray, k: int):
    '''
    Given an expected motif sort configurations into n groups to maximise convergence. 
    '''
    
    print('Beginning to sort configs for smart convergence. This process can be costly')

    group_size = len(configs)/ngroups
    indicies = np.zeros((ngroups,int(len(configs)/ngroups)),dtype=np.int16)
    counts = np.zeros(ngroups,dtype=np.int16)
    available_mask = np.ones(len(configs), dtype=bool) #mask used configs.
    all_pdds = [PDD(c.positions, c.cell, k) for c in configs]
    # first find the starting point for each group.
    init_emds = np.zeros((len(configs)))
    for i,config in enumerate(configs):
        motif_config = return_motif_config(config,expected_motif)
        PDD1 = all_pdds[i]
        PDD2 = PDD(motif_config.positions,motif_config.cell,k)
        init_emds[i] = EMD(PDD1,PDD2) 
    
    end_points = np.argpartition(init_emds, ngroups-1)[:ngroups]
    available_mask[end_points] = False
    seen_configs = [configs[i] for i in end_points]
    indicies[:,0] = end_points
    
    # iterative expansion of groups
    while np.any(counts < group_size-1):
        progress = np.sum(counts+1)/len(configs) * 100
        print(f"\rProgress: {(round(progress,1))}%", end="")

        config_to_append = None
        config_to_push_back = None
        min_emd = 1 # max value of the emd.
        for i,config in enumerate(configs):
            if not available_mask[i]:
                continue
            for j,end_index in enumerate(end_points):
                if counts[j] != (group_size-1): 
                    PDD1 = all_pdds[i]
                    PDD2 = all_pdds[end_index]
                    emd = EMD(PDD1,PDD2)
                    if emd <= min_emd:
                        min_emd = emd
                        config_to_append = i
                        config_to_push_back = j
        
        # update
        end_points[config_to_push_back] = config_to_append
        available_mask[config_to_append] = False
        seen_configs.append(configs[config_to_append])     
        counts[config_to_push_back] += 1
        indicies[config_to_push_back][counts[config_to_push_back]] = config_to_append
        
    print('Done grouping')
    
    return indicies
    

def return_motif_config(config: Atoms, motif: np.ndarray):
    '''
    Some positions may be wrapped to larger cell sizes. 
    '''
    motif_config = config.copy()
    lattice_vectors = np.array(config.cell)
    #first search surrounding space in case config positions are wrapped. (may need to be edited for non-square cells?)
    motif_cart_extended = []
    for i,j,k in product(range(0,2),range(0,2),range(0,2)):
        for motif_pos in motif:
            pos = (motif_pos + np.array([i,j,k]))
            pos_cart = pos[0] * lattice_vectors[0] + pos[1] * lattice_vectors[1] + pos[2] * lattice_vectors[2]
            motif_cart_extended.append(pos_cart)
    # notice this is very similar to that used in mlipts.codes.vasp.set_magmom, could be generalized. 
    A = config.positions
    B = np.array(motif_cart_extended)
    diff = A[:, None, :] - B[None, :, :]  
    dist2 = np.sum(diff**2, axis=2)       
    closest_indices = np.argmin(dist2, axis=1)
    motif_cart = B[closest_indices]
    motif_config.set_positions(motif_cart)

    return motif_config
    