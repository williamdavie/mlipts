'''

@author William Davie

Implementation of the Pointwise Distance Distribution for similarity assessment of material data sets.

Following the work of Daniel Widdowson and Vitaliy Kurlin (2021): https://arxiv.org/abs/2108.04798v3
'''

import numpy as np
from itertools import product
from typing import Iterator
from scipy.spatial import KDTree
from ase.io import read,write 
from ase import Atoms


def atoms_configs_PDDs(configs: list[Atoms], k: int):
    '''
    Computes the pointwise distance distribution for all entries in a list of ase.Atoms objects
    '''
    PDDs = []
    
    for i in range(len(configs)):
        cell = np.array(configs[i].cell)
        motif = np.array(configs[i].positions)
        PDDs.append(PDD(motif,cell,k))
    return PDDs
    

def PDD(motif: np.ndarray,
        cell: np.ndarray,
        k: int) -> np.ndarray: 
    '''
    Computes the pointwise distance distribution, implemented according to the pseudocode of (Widdowson, et al. 2021)
    
    Parameters
    ----------
    motif : :class:`numpy.ndarray` 
        Nx3 dimensional motif describing the atomic positions in a periodic lattice.
    cell : :class:`numpy.ndarray` 
        lattice vectors.
    k : int
        number of neighbours to consider
    
    Returns
    -------
    PDD: :class:`numpy.ndarray` 
        an array of shape (atoms in motif,k)
    
    References
    ----------
    [1] Widdowson, Daniel, and Vitaliy Kurlin. "Pointwise distance distributions for detecting near-duplicates in large materials databases." arXiv preprint arXiv:2108.04798 (2021).
    '''
    
    S = np.empty((0,3)) # initially S = M
    g = point_generator(motif,cell,k)
    
    # minimum required points
    while len(S) < k+1:
        points = next(g)
        S = np.append(S,points,axis=0)
    tree = KDTree(S)
    #init D(S,M;k)
    D_, inds = tree.query(motif, k+1)
    
    while True:
        
        D_prev = np.copy(D_)
        points = next(g)
        S = np.append(S,points,axis=0)
        
        tree = KDTree(S)
        D_, inds = tree.query(motif,k+1)
        
        if np.array_equal(D_,D_prev):
            break
        
    PDD_collapsed = collapse(D_)
    PDD_final = PDD_collapsed[np.lexsort(PDD_collapsed.T[::-1])]
    
    return PDD_collapsed

def point_generator(motif: np.ndarray,
              cell: np.ndarray,
              k: int) -> Iterator[np.ndarray]:
    '''
    An iterator to sequentially add points to set of atomic positions S. Looping through the motif and lattice vectors and sequentially adding atoms.
    '''
    
    for upper_bound in range(k+1): # move 1 shell up
        # will never need to extend the lattice more than ±(ka + kb + kc).
        points = np.empty((0,3),dtype=np.float32)
        for (i,j,k) in product(range(-upper_bound,upper_bound+1),repeat=3):
            if max(abs(i),abs(j),abs(k)) != upper_bound:
                continue
    
            new_points = motif + (cell[0] * i) + (cell[1] * j) + (cell[2] * k)
            points = np.append(points,new_points,axis=0)

        yield points
        
        
def collapse(D: np.ndarray,tol: int=3) -> np.ndarray:
    '''
    Collapses repeated entries in the matrix D(S,M;k). 
    Weights are associated with repeated entries (l/m, where l is the number of repetitions)
    '''
    
    nAtoms,k_p1 = D.shape

    D_rounded = np.round(D, decimals=tol) # tolerance level for equality. 
    D_unique,counts = np.unique(D_rounded,axis=0,return_counts=True)
    multipliers = np.ones_like(counts,dtype=np.float32)
    
    for i,count in enumerate(counts):
        multipliers[i] = count/nAtoms
            
    D_unique[:,0] = multipliers
    
    
    return D_unique