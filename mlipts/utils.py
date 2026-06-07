"""
MLIPTS utilities.
"""

from itertools import product

import ase
import ase.build
import ase.neighborlist
import numpy as np
import scipy.spatial

# ------------------------------Reference positions------------------------------------


def get_scaled_reference_positions(config: ase.Atoms, equilibrium_config: ase.Atoms):
    """
    Returns the scaled positions of a input config reformatted to align with the equilibrium configuration (high-symmetry lattice structure)
    """

    cell = config.cell
    minimal_cell = equilibrium_config.cell

    scaling_factor = int(
        np.linalg.det(cell) / np.linalg.det(minimal_cell)
    )  # intergers 1,2,3 .. etc

    # should be close to 1
    volume_scaling_factor = (np.linalg.det(cell) / scaling_factor) / np.linalg.det(
        minimal_cell
    )
    expanding_vectors = (volume_scaling_factor) ** (1 / 3) * minimal_cell

    # reference position, this is jubious but the main thing is that we select the correct element as reference
    reference_element = get_reference_element(equilibrium_config)
    reference_positions = get_reference_positions(config)
    reference_position = reference_positions[reference_element]

    input_pos_ref = config.get_positions() - reference_position
    input_pos_scaled = np.zeros_like(input_pos_ref)

    for i, pos in enumerate(input_pos_ref):
        input_pos_scaled[i] = np.linalg.solve(expanding_vectors, pos)

    return input_pos_scaled


def get_reference_positions(atoms: ase.Atoms):

    reference_positions = {}
    dict_keys = []
    atoms_symbols = atoms.get_chemical_symbols()
    atoms_positions = atoms.get_positions()

    for element in atoms_symbols:
        if element not in dict_keys:
            dict_keys.append(element)

            elem_positions = atoms_positions[[s == element for s in atoms_symbols]]
            closest_pos = elem_positions[
                np.argmin(np.linalg.norm(elem_positions, axis=1))
            ]

            reference_positions[str(element)] = closest_pos

    return reference_positions


def get_reference_element(atoms: ase.Atoms):

    pos = atoms.positions
    idx = np.argmin(np.linalg.norm(pos, axis=1))

    return atoms[idx].symbol


# -------------------------------------------------------------------------------------


def sort_configs_by_volume(configs: list[ase.Atoms]) -> list[ase.Atoms]:
    """
    Sort a set of configurations by volume
    """

    volumes = np.array([c.get_volume() for c in configs])
    indicies = np.argsort(volumes)
    return [configs[i] for i in indicies]


# -------------------------------Supercells---------------------------------------------


def match_config_to_dir(config: ase.Atoms, supercell_dict: dict) -> str:
    """
    Given a directory, labelled by a set of lattice vectors, match the config to a directory.
    """

    dirs = list(supercell_dict.keys())
    values = list(supercell_dict.values())
    min_val = 100
    index = 0

    for i, value in enumerate(values):

        dif = abs(np.linalg.norm(value - np.array(config.cell)))

        if dif < min_val:
            min_val = dif
            index = i

    if min_val > np.min(values[index]) / 2:
        print("Warning <!>: Did not find a matching supercell in the input dictionary")

    return dirs[index]


def fetch_supercell_motif(motif: np.ndarray, supercell_dims: np.ndarray):

    Nx, Ny, Nz = supercell_dims

    supercell_motif = []
    for i, j, k in product(range(0, Nx), range(0, Ny), range(0, Nz)):
        for pos in motif:
            new_pos = pos + np.array([i, j, k])
            supercell_motif.append(new_pos / supercell_dims)

    return np.array(supercell_motif)


# ------------------------------Retrieving supercell matricies--------------------------------------


def get_supercell_matrix(cell, minimal_cell):
    """
    Get the transformation between the original high symmetry primitive lattice vectors and the input cell.
    """

    scaling_factor = round(
        np.linalg.det(cell) / np.linalg.det(minimal_cell)
    )  # intergers 1,2,3 .. etc

    # should be close to 1
    volume_scaling_factor = (np.linalg.det(cell) / scaling_factor) / np.linalg.det(
        minimal_cell
    )

    expanded_minimal_cell = (volume_scaling_factor) ** (1 / 3) * minimal_cell

    # H' = SH, solve for S
    S = np.linalg.solve(expanded_minimal_cell.T, cell.T).T

    return S


def retrieve_standard_supercell(cell, minimal_cell, tol):
    """
    Given a rotated non-diagonal supercell and the primitive lattice vectors. Calculate the interger supercell matrix.

    This is achieved by mapping the angles and norm ratios between vectors for all possible supercell matricies with the expected determinant.
    """

    supercell_det = round(np.linalg.det(cell) / np.linalg.det(minimal_cell))

    cell_parameters = get_magnitudes_and_angles(cell)

    possible_supercells = get_possible_supercells(supercell_det)

    for supercell in possible_supercells:

        supercell_parameters = get_magnitudes_and_angles(supercell)

        if np.allclose(supercell_parameters, cell_parameters, atol=tol):

            parameter = np.linalg.norm(cell) / np.linalg.norm(supercell)

            result = parameter * supercell

            return result
    return None


def get_possible_supercells(det: int):
    """
    given the desired matrix determinant, return all possible upper triangular HNF supercell matricies
    """

    # fetch all diagonal elements
    diagonal_elements = []
    for i in range(1, det + 1):
        if det % i == 0:
            remaining = det // i
            for j in range(1, remaining + 1):
                if remaining % j == 0:
                    k = remaining // j
                    diagonal_elements.append((i, j, k))

    # define all given a diagonal element
    supercells = []
    for diag in diagonal_elements:
        S00, S11, S22 = diag
        for n in range(0, diag[1]):
            S01 = n
            for p, q in product(range(0, diag[2]), range(0, diag[2])):
                S02, S12 = p, q

                S = np.array([[S00, S01, S02], [0, S11, S12], [0, 0, S22]])

                supercells.append(S)

    return supercells


def get_magnitudes_and_angles(cell):

    alpha = vectors_angle(cell[1], cell[2])
    beta = vectors_angle(cell[0], cell[2])
    gamma = vectors_angle(cell[0], cell[1])

    mag_array = np.array(
        [np.linalg.norm(cell[0]), np.linalg.norm(cell[1]), np.linalg.norm(cell[2])]
    )
    mag_array /= np.min(mag_array)

    return np.array([alpha, beta, gamma, mag_array[0], mag_array[1], mag_array[2]])


def vectors_angle(vec1, vec2):

    cosAngle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))

    return np.arccos(cosAngle)


def return_motif_config(
    config: ase.Atoms, equilibrium_config: ase.Atoms, atol: float = 1e-4
) -> ase.Atoms:
    """
    Some positions may be wrapped to larger cell sizes.
    """
    motif_config = config.copy()

    cell = config.cell
    minimal_cell = equilibrium_config.cell
    # create a high symmetry configuration for an expanded cell
    expanded_config = equilibrium_config.copy()

    scaling_factor = round(np.linalg.det(cell) / np.linalg.det(minimal_cell))
    volume_scaling_factor = (np.linalg.det(cell) / scaling_factor) / np.linalg.det(
        minimal_cell
    )
    expanded_cell = (volume_scaling_factor) ** (1 / 3) * minimal_cell

    expanded_config.set_cell(expanded_cell)
    expanded_config.set_scaled_positions(equilibrium_config.get_scaled_positions())

    # now expand this cell to the input supercell.

    S = get_supercell_matrix(cell, expanded_cell)

    if np.allclose(S, np.round(S), atol=atol):
        supercell_motif = ase.build.make_supercell(expanded_config, np.round(S))
    else:
        raise ValueError(
            "cannot set magmom for configurations with non-interger supercell matricies, maybe need to use a conversion before using as input data."
        )

    positions_extended = []
    for i, j, k in product(range(-1, 2), range(-1, 2), range(-1, 2)):
        for pos in supercell_motif.get_positions():
            new_pos = pos + np.array([i, j, k]) @ supercell_motif.cell
            positions_extended.append(new_pos)

    # notice this is very similar to that used in mlipts.codes.vasp.set_magmom, could be generalized.
    A = config.get_positions()
    B = np.array(positions_extended)
    diff = B[None, :, :] - A[:, None, :]
    dist2 = np.sum(diff**2, axis=2)
    closest_indices = np.argmin(dist2, axis=1)
    positions_new = B[closest_indices]

    motif_config.set_positions(positions_new)

    return motif_config


def map_properties(
    configA: ase.Atoms, configB: ase.Atoms, supercell_atol: float = 1e-4
):
    """
    maps the features of configB to configA independant of supercell.

    Only copies across features that aren't present in configB.

    Again we have lots of parallels between this function, return_motif_config and codes.vasp.set_magmom).
    """
    # configB is the expected minimal cell, it may be the same size as B.

    configA.wrap()

    minimal_cell = np.array(configB.cell)

    cell = np.array(configA.cell)

    try:
        configB_scaled_pos = configB.get_scaled_positions()
    except Exception as exc:
        raise ValueError("config A must contain positions") from exc

    scaling_factor = round(np.linalg.det(cell) / np.linalg.det(minimal_cell))
    volume_scaling_factor = (np.linalg.det(cell) / scaling_factor) / np.linalg.det(
        minimal_cell
    )
    expanded_cell = (volume_scaling_factor) ** (1 / 3) * minimal_cell

    # scale up the positions, all other properties i.e. magnetic moment remain the same.
    expanded_configB = configB.copy()
    expanded_configB.set_cell(expanded_cell)
    expanded_configB.set_scaled_positions(configB_scaled_pos)

    # now expand this cell to the input supercell.
    S = get_supercell_matrix(cell, expanded_cell)

    if np.allclose(S, np.round(S), atol=supercell_atol):
        supercell_configB = ase.build.make_supercell(expanded_configB, np.round(S))
    else:
        raise ValueError(
            "cannot copy data across configurations with non-interger supercell matricies, maybe need to use a conversion before using as input data."
        )

    supercell_configB.wrap()

    # retrieve desired arrays:
    final_arrays = {}
    for Bkey, B_property in supercell_configB.arrays.items():
        if Bkey not in configA.arrays.keys():
            if isinstance(B_property, np.ndarray):
                shape = (len(configA),) + B_property.shape[1:]
                final_arrays[Bkey] = np.zeros(shape)
            else:
                final_arrays[Bkey] = np.zeros(len(configA))

    # Now match arrays according to minimum distance and element.
    A = configA.positions
    B = configB.positions
    all_diff = B[None, :, :] - A[:, None, :]
    all_dist2 = np.sum(all_diff**2, axis=2)

    for i, atom in enumerate(configA):
        symbol = atom.symbol
        element_indices = [
            index for index, val in enumerate(configB.symbols) if val == symbol
        ]
        if not element_indices:
            raise ValueError("reference configuration has missing elements. ")

        B_this_element = B[np.array(element_indices), :]
        diff = B_this_element - A[i]
        dist2 = np.sum(diff**2, axis=1)
        closest_index = np.argmin(dist2)
        true_index = np.where(dist2[closest_index] == all_dist2)[1][0]
        for key, final_property in final_arrays.items():
            final_property[i] = supercell_configB.arrays[key][true_index]

    for key, final_property in final_arrays.items():
        configA.arrays[key] = final_property

    return configA


# -------------------------------Defects---------------------------------------------


def generate_defect(
    config: ase.Atoms,
    targets: dict[str, int],
    defect_type: str = "schottky",
    order: int = 1,
) -> ase.Atoms:
    """
    Generates a config that contains a defect.

    Parameters
    ----------
    order: int
        for defect_type 'schottky' this gives nth order neighest neighbours of the highest charged element.
        for defect_type 'frenkel' this gives the promixity of the involved atom and it's original initistitual pair.
    """
    config.wrap()
    if defect_type == "schottky":
        config = generate_schottky_defect(config, targets, neighbour_order=order)
    elif defect_type == "frenkel":
        config = generate_frenkel_defect(config, targets, proximity=order)

    return config


def generate_schottky_defect(
    config: ase.Atoms, targets: dict[str, int], neighbour_order: int
):
    """generate schottky defect

    Parameters
    ----------
    config:
        atomic configuration to generate defect in
    targets: dict[str, int]
        target elements to remove. For example for a 4+ cation (X) and 2- anion (A) pair the dict for a schottky defect would be {'X' : 1, 'A' : 2}.
    neighbour_order: int
        the number of unit cells seperating the highest charge element and the rest of the defect.

    Returns
    -------
    config:
        atoms with defect removed.
    """

    # find smallest atom count (should have highest charge).
    minimum_count = min(targets.values())
    minimum_key = [key for key in targets if targets[key] == minimum_count][
        0
    ]  # if more than one value just select the first.

    # choose largest elements
    old_config = config.copy()
    symbols = np.array(config.get_chemical_symbols(), dtype="U10")
    removed_largest = []
    removed = []
    min_key_indices = np.where(symbols == minimum_key)[0]
    for idx in min_key_indices[:minimum_count]:
        removed_largest.append(idx)
        removed.append(idx)

    # remove nearest neighbours of highest charge element.
    surrounding_elements = [key for key in targets if targets[key] != minimum_count]
    for j in surrounding_elements:
        to_remove = find_nearest_neighbours(
            old_config,
            removed_largest,
            Nneighbours=targets[j],
            target_element=j,
            order=neighbour_order,
        )
        removed.extend(to_remove)

    for idx in sorted(removed, reverse=True):
        del config[idx]

    return config


def generate_frenkel_defect(config: ase.Atoms, targets: dict[str, int], proximity: int):
    """generate frenkel defect"""

    symbols = np.array(config.get_chemical_symbols(), dtype="U10")

    for target_element in targets.keys():
        atom_indices = np.where(symbols == target_element)[0]
        count = targets[target_element]

        # this sets a preference for the 'largest availible' interstitial.
        interstitial_sites = find_interstitial_sites(
            config,
            count + (count * proximity),
        )

        for i in range(count):
            # finds the (n+1)th furtherest interstitual site using proximity parameter
            orig_pos = config.positions[atom_indices[i]].copy()
            distances = scipy.spatial.distance.cdist([orig_pos], interstitial_sites)[0]
            sorted_indices = np.argsort(distances)
            print(f"chosen site: {interstitial_sites[sorted_indices][-1]}")
            config.positions[atom_indices[i]] = interstitial_sites[sorted_indices][-1]

    return config


def find_nearest_neighbours(
    config: ase.Atoms,
    central_atoms: list[int],
    Nneighbours: int,
    target_element: str = "O",
    order: int = 1,  # 1 - nearest neighbour, 2 - next neighest neighbour etc
):
    """Finds nearest N (nth order) neighbours of central_atoms (indicies of the config.)"""

    i, j, d = ase.neighborlist.neighbor_list("ijd", config, cutoff=np.max(config.cell))

    symbols = np.array(config.get_chemical_symbols(), dtype="U10")
    is_target = symbols[j] == target_element

    for atom in central_atoms:
        is_central = i == atom
        mask = is_central & is_target
        neighbor_indices = j[mask]
        neighbor_distances = d[mask]

        if len(neighbor_distances) == 0:
            raise ValueError("No neighbours of this target  element found")

        ordered_idx = np.argsort(neighbor_distances)
        sorted_indices = neighbor_indices[ordered_idx]
        sorted_dists = neighbor_distances[ordered_idx]

        norm_dists = np.round(sorted_dists / np.min(sorted_dists), 1)
        unique_norm = np.unique(norm_dists)

        if order <= len(unique_norm):
            shell_mask = norm_dists == unique_norm[order - 1]
            result_indices = sorted_indices[shell_mask]
        else:
            raise ValueError(
                f"Cannot find the {order}th nearest neighbors in this supercell."
            )

    return result_indices[0:Nneighbours]


def insert_cluster_to_interstitual(config: ase.Atoms, cluster: ase.Atoms):
    """Centres a cluster at an intersititial_site"""

    interstitial_site = find_interstitial_sites(config, 1, len(config))

    # find cluster centre,
    cluster_positions = cluster.get_positions()
    centroid = np.mean(cluster_positions, axis=0)
    centered_positions = cluster_positions - centroid
    # normalize positions around that centre
    centered_positions = cluster_positions - centroid
    # place the cluster at initerstitial site
    new_positions = centered_positions + interstitial_site
    for i, atom in enumerate(cluster):
        new_atom = ase.Atom(symbol=atom.symbol, position=new_positions[i])
        config.append(new_atom)

    return config


def insert_element_to_interstitial(config: ase.Atoms, element: str, count: int = 1):
    """Adds N elements to N interstitial sites"""

    interstitial_sites = find_interstitial_sites(config, count)

    for site in interstitial_sites:
        new_atom = ase.Atom(element, position=site)
        config.append(new_atom)

    return config


def find_interstitial_sites(atoms: ase.Atoms, count: int, grid_density: int = 50):
    """Basic function to find an interstitial site."""

    pos = atoms.get_positions()
    cell = atoms.get_cell()

    x = np.linspace(0, 1, grid_density, endpoint=False)
    grid_points_frac = np.array(np.meshgrid(x, x, x)).T.reshape(-1, 3)
    grid_points_cart = np.dot(grid_points_frac, cell)

    shifts = np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])).T.reshape(-1, 3)
    cart_shifts = np.dot(shifts, cell)
    expanded_pos = np.vstack([pos + shift for shift in cart_shifts])

    tree = scipy.spatial.KDTree(expanded_pos)
    distances, _ = tree.query(grid_points_cart)

    best_indices = np.argsort(distances)[::-1]
    best_sites = grid_points_cart[best_indices[0:count]]

    print(f"Found interstitual sites: {best_sites}")

    return best_sites
