"""
utils for analysing a constructed dataset.
"""

import ase
import ase.build
import numpy as np
from ase.geometry import get_distances
from prettytable import PrettyTable
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

from mlipts.utils import get_supercell_matrix

# -------------- Analyse a training data set ----------------------


def identify_neutral_defect(
    config: ase.Atoms,
    equalibrium_config: ase.Atoms,
    target: dict[str, int],
    defect_type: str,
    frenkel_tol_multiplier: float = 0.3,
) -> bool:
    """returns True of defect with selected targets and given defect type is found in Atoms"""

    supercell_matrix = np.round(
        get_supercell_matrix(config.cell, equalibrium_config.cell), 0
    )
    equalibrium_config = ase.build.make_supercell(equalibrium_config, supercell_matrix)
    equalibrium_config.set_cell(config.cell, scale_atoms=True)
    lattice_param = (config.get_volume() / np.linalg.det(supercell_matrix)) ** (1 / 3)

    found = False
    if defect_type == "frenkel":
        found = identify_frenkel(
            config, equalibrium_config, target, lattice_param, frenkel_tol_multiplier
        )
    if defect_type == "schottky":
        found = identify_schottky(config, equalibrium_config, target)

    return found


def identify_schottky(
    config: ase.Atoms, equalibrium_config: ase.Atoms, target: dict[str, int]
):

    # easy to identify - any additional features from bulk should be removed.
    if len(equalibrium_config) - sum(target.values()) == len(config):
        return True
    else:
        return False


import ase
import numpy as np
from scipy.spatial import cKDTree


def get_non_equalibrium_atoms(
    config: ase.Atoms, equilibrium_config: ase.Atoms, cutoff: float
):
    """
    Identifies interstitials using Minimum Image Convention (MIC)
    compatible with non-orthogonal (triclinic) cells.
    """
    interstitials = []

    config.wrap()
    equilibrium_config.wrap()

    symbols = set(config.get_chemical_symbols()) | set(
        equilibrium_config.get_chemical_symbols()
    )

    for symbol in symbols:
        c_sub = config[config.symbols == symbol]
        e_sub = equilibrium_config[equilibrium_config.symbols == symbol]

        if len(c_sub) == 0 or len(e_sub) == 0:
            continue
        c_pos = c_sub.get_positions()
        e_pos = e_sub.get_positions()

        from ase.geometry import get_distances

        _, dists = get_distances(c_pos, e_pos, cell=config.get_cell(), pbc=True)

        min_dists = np.min(dists, axis=1)
        for i, d in enumerate(min_dists):
            if d > cutoff:
                interstitials.append(c_sub[i])

    return interstitials


def identify_frenkel(
    config: ase.Atoms,
    equilibrium_config: ase.Atoms,
    target: dict[str, int],
    lattice_param: float,
    frenkel_tol_multiplier: float = 0.4,
):

    config.wrap()
    equilibrium_config.wrap()

    if len(config) != len(equilibrium_config):
        return False

    cutoff = lattice_param * frenkel_tol_multiplier
    interstitials = get_non_equalibrium_atoms(config, equilibrium_config, cutoff)

    # In a fixed-count system, if an atom is > cutoff from home, it means another site is empty (a vacancy).
    for symbol, count in target.items():
        count_found = sum(1 for a in interstitials if a.symbol == symbol)

        if count_found < count:
            return False

    return True


def configuration_distribution_table(atoms: list[ase.Atoms]):
    """Given a set of atoms returns a table showing total number of configs for each number of atoms
    i.e.:

    N atoms : N configs
    24 : 1000
    48 : 500
    """
    Natoms = np.array([atom.get_number_of_atoms() for atom in atoms], dtype=np.int64)

    unique_Natoms, totals = np.unique(Natoms, return_counts=True)

    table = PrettyTable()

    table.field_names = ["N atoms", "Total configurations"]

    for i, count in enumerate(unique_Natoms):
        table.add_row([f"{count}", f"{totals[i]}"])

    return table


def supercell_distribution_table(atoms: list[ase.Atoms], minimal_atoms: ase.Atoms):
    """Given a set of atoms returns a table showing total number of configs for each supercell"""

    supercells = np.array(
        [
            np.round(get_supercell_matrix(atom.cell, minimal_atoms.cell)).astype(
                np.int32
            )
            for atom in atoms
        ]
    )

    unique_flat = np.unique(supercells.reshape(supercells.shape[0], -1), axis=0)
    unique_supercells = unique_flat.reshape(-1, 3, 3)

    _, counts = np.unique(
        supercells.reshape(supercells.shape[0], -1), axis=0, return_counts=True
    )

    table = PrettyTable()

    table.field_names = ["Super cell", "Total configurations", "Det"]

    for i, supercell in enumerate(unique_supercells):
        table.add_row(
            [
                f"{tuple(map(tuple, supercell))}",
                f"{counts[i]}",
                f"{np.linalg.det(supercell)}",
            ]
        )

    return table
