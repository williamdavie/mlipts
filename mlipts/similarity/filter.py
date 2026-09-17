"""
@author William Davie

filter functions
"""

import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from scipy.cluster import hierarchy

from mlipts.similarity.emd import EMD, cached_EMD
from mlipts.similarity.pdd import atoms_configs_PDDs


def filter_by_emd(
    configs: list[Atoms],
    tol: float,
    k: int,
    show_dendrograms: bool = False,
    show_results: bool = False,
):
    """
    filters a set of configurations by earth movers distance according to a tolerance
    """

    print(
        "---------------------------------------------------------------------------------------"
    )
    print(
        f"Beginning to filter {len(configs)} configurations by EMD, this process can be costly."
    )
    print(
        "---------------------------------------------------------------------------------------"
    )

    PDDs = atoms_configs_PDDs(configs, k)
    emds = []
    emd_cache = {}
    emds_filtered = []

    remove_ = set()
    print("caching implemented")
    for i in range(len(PDDs)):
        print(f"\rProgress: {int(100*i/len(PDDs))}%", end="")
        if i in remove_:
            continue
        for j in range(i + 1, len(PDDs)):
            emd = cached_EMD(i, j, PDDs, emd_cache)
            if show_results:
                print(emd)
            if j in remove_:
                continue
            if emd <= tol:
                remove_.add(j)
            else:
                emds_filtered.append(emd)
    print("\n")

    inds = [i for i in range(len(PDDs)) if i not in remove_]

    return [configs[i] for i in inds], inds
