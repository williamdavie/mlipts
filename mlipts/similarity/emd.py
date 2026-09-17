"""

@author: William Davie

Computes the Earth Movers distance to access similarity between two PDDs.

Following the work of Daniel Widdowson and Vitaliy Kurlin (2021): https://arxiv.org/abs/2108.04798v3.

Calculation of earth movers distance is an optimization problem, see Sec.4 of reference, there are a number of algorithms to solve it.
This implementation takes advantage of scipy.optimize.linprog noting the oppitunity to add the network simplex problem specific algorithm.
"""

import matplotlib.pyplot as plt
import numpy as np
from ase.io import read, write
from scipy.cluster import hierarchy
from scipy.optimize import linprog
from scipy.spatial.distance import cdist


def EMD(pdd_S: np.ndarray, pdd_Q: np.ndarray) -> float:
    """
    Calculates the earth movers distance given two PDDs (Pointwise pair distributions) on metric spaces S and Q.

    Parameters
    ----------
    pdd_S : :class:`numpy.ndarray`
        The PDD for metric space S

    pdd_Q  : :class:`numpy.ndarray`
        The PDD for metric space Q

    Returns
    -------
    emd : float
        The Earth Movers Distance EMD to describe similarity between two distributions.

    """
    # separate pdd into weights and d_k values
    wS, wQ, dist_S, dist_Q = pdd_S[:, 0], pdd_Q[:, 0], pdd_S[:, 1:], pdd_Q[:, 1:]
    # compute a cost function (distance metric between PDDs)
    D_SQ = cdist(dist_S, dist_Q, metric="chebyshev")

    # calculate emd via scipy - leaving oppitunity to add other algorithms as functions.
    return emd_v_scipy(wS, wQ, D_SQ)


def cached_EMD(i: int, j: int, pdds: list[np.ndarray], emd_cache: dict) -> float:
    """
    Caching the emd (trading memory for speed up), stores the emd in a dictionary to be re-used when required.
    """

    key = (i, j) if i <= j else (j, i)
    if key not in emd_cache:
        emd_cache[key] = EMD(pdds[i], pdds[j])
    else:
        pass
    return emd_cache[key]


def EMD_hierarchy(PDDs: np.ndarray) -> None:
    """
    Calculates the Earth Movers Distance across a set of PDDs and plots a corropsonding 'similarity' hierarchy.

    Copyright (C) 2025 Daniel Widdowson
    """

    emds = []
    for i in range(len(PDDs)):
        print(f"\rProgress: {int(100*i/len(PDDs))}%", end="")
        for j in range(i + 1, len(PDDs)):
            emd = EMD(PDDs[i], PDDs[j])
            emds.append(emd)
    emds = np.array(emds)
    print(len(emds))

    Z = hierarchy.linkage(emds)
    dn = hierarchy.dendrogram(Z)
    plt.show()

    return None


def emd_v_scipy(wS: np.ndarray, wQ: np.ndarray, D_SQ: np.ndarray) -> float:
    """
    Uses scipy.optimize.linprog to calculate the earths movers distance.
    wS and wQ are weights to define constraints on f_ij.
    D_SQ is the distance (cost function).
    """
    mS, mQ = D_SQ.shape
    n_vars = mS * mQ
    n_constraints = mS + mQ

    A_eq = np.zeros(shape=(n_constraints, n_vars))

    for i in range(mS):
        A_eq[i, i * mQ : (i + 1) * mQ] = 1

    for j in range(mQ):
        A_eq[mS + j, j:n_vars:mQ] = 1

    b_eq = np.concatenate([wS, wQ])

    return linprog(
        D_SQ.flatten(), A_eq=A_eq, b_eq=b_eq, bounds=[(0, 1)], method="highs"
    ).fun


def compute_emd_pair(args):
    i, j, pdd1, pdd2 = args
    return i, j, EMD(pdd1, pdd2)
