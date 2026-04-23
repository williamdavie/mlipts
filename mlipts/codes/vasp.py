"""
@author: William Davie

File containing vasp specific functionality. Used to build many vasp calculations.
"""

import shutil
import sys
from itertools import product
from pathlib import Path

import ase
import ase.build
import ase.io
import numpy as np

from mlipts import utils


def build_vasp_calculation(
    vasp_base_dir: str, config: ase.Atoms, calc_name: str, outdir: str
) -> str:
    """
    Builds a vasp calculation directory for a given atomic configuration.

    Parameters
    ----------
    vasp_base_dir: str
        path to base directory, should contain POTCAR, KPOINTS and INCAR.
    config: :class:`ase.Atoms`
        atomic configuration.
    outname: stls
        name of calculation directory
    outdir: str
        output path of calculation directory.

    Returns
    -------
    new_calc_dir: str
        vasp directory generated in outdir.
    """

    new_calc_dir = outdir + "/" + f"{calc_name}"
    shutil.copytree(vasp_base_dir, new_calc_dir, dirs_exist_ok=True)
    ase.io.write(
        new_calc_dir + "/POSCAR", config, format="vasp", vasp5=True, direct=True
    )

    # set magnetic moment
    if "magnetic_moment" in config.arrays.keys():
        print("Found magnetic moment in input file, now setting MAGMOM")
        magmom_str = "MAGMOM = "
        for momoment in config.arrays["magnetic_moment"]:
            magmom_str += f"{momoment[0]} {momoment[1]} {momoment[2]} "

        writeMAGMOM(f"{new_calc_dir}/INCAR", new_magmom_str=magmom_str)

    if "Uvalue" in config.arrays.keys():
        print("Found Uvalue in input file, now setting LDAUU")
        symbols = config.get_chemical_symbols()
        map_ = dict.fromkeys(symbols)

        for i, value in enumerate(config.arrays["Uvalue"]):
            map_[symbols[i]] = value

        Uvalues = [value for key, value in map_.items()]

        set_Uvalue(new_calc_dir, Uvalues)

    return new_calc_dir


def write_POSCAR_str(config: ase.Atoms) -> str:
    """
    writes a POSCAR string given an atomic configuration.
    """

    poscar = "System\n 1.0\n"

    cell = np.array(config.cell)
    poscar += f" {cell[0,0]} {cell[0,1]} {cell[0,2]}\n {cell[1,0]} {cell[1,1]} {cell[1,2]}\n {cell[2,0]} {cell[2,1]} {cell[2,2]}\n"

    type_list = list(config.symbols)

    # set can be unordered so can't use set(config.symbols)
    type_labels = []
    for i in type_list:
        if i not in type_labels:
            type_labels.append(i)  # only way i can see to gaurentee order?
    for element in type_labels:
        poscar += f" {element} "
    poscar += "\n"

    for element_type in type_labels:
        count = config.symbols.count(element_type)
        poscar += f" {count} "
    poscar += "\nCartesian\n"

    for pos in config.get_positions():
        poscar += f"{pos[0]} {pos[1]} {pos[2]}\n"

    return poscar


def append_vasp_calc_to_database(
    database_file: str,
    vasp_dir: str,
    save_magmoms: bool = True,
    save_charge: bool = True,
    store_failed: bool = True,
    store_Uvalue: bool = True,
) -> None:

    # read atom data.
    atoms = ase.io.read(f"{vasp_dir}/vasprun.xml")

    with open(f"{vasp_dir}/OUTCAR", "r", encoding="utf-8") as f:
        outcar_str = f.read()
    if "aborting loop EDIFF was not reached (unconverged)" in outcar_str:
        print(f"Self consistency failed in {vasp_dir}, not saving data.")
        if store_failed:
            ase.io.write("failed_to_converge.xyz", atoms, format="extxyz", append=True)
        return None

    if save_magmoms:
        magmoms = fetch_OUTCAR_magnetization(f"{vasp_dir}/OUTCAR", len(atoms))
        atoms.set_array("magnetic_moments", magmoms)
    if save_charge:
        charges = fetch_OUTCAR_total_charge(f"{vasp_dir}/OUTCAR", len(atoms))
        atoms.set_array("total_charge", charges)
    if store_Uvalue:
        u_vals = fetch_OUTCAR_Uvalues(f"{vasp_dir}/OUTCAR", atoms)
        atoms.set_array("Uvalue", u_vals)

    ase.io.write(database_file, atoms, format="extxyz", append=True)
    return None


def fetch_OUTCAR_total_charge(outcar: str, atom_count: int):
    """Reads total charge from OUTCAR"""
    return read_OUTCAR_atom_totals(outcar, "total charge", atom_count)


def fetch_OUTCAR_magnetization(outcar: str, atom_count: int):
    """Reads magnetic_moments from OUTCAR"""

    drcts = ["x", "y", "z"]

    magmoms = np.zeros((atom_count, 3))

    for i, drct in enumerate(drcts):
        magmoms[:, i] = read_OUTCAR_atom_totals(
            outcar, f"magnetization ({drct})", atom_count
        )

    return magmoms


def read_OUTCAR_atom_totals(outcar: str, label: str, atoms_count: int):
    """Reads atom data from outcar given a label"""

    with open(outcar, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # want to store the last instance of this data.
    location = 0
    for i, line in enumerate(reversed(lines)):
        if label in line:
            location = len(lines) - i
            break

    data = []
    for line in lines[location + 3 : location + 3 + atoms_count]:
        vals = line.split()
        data.append(float(vals[-1]))

    return np.array(data)


def fetch_OUTCAR_Uvalues(outcar: str, atoms: ase.Atoms):

    # map chemical symbols to U value
    map_ = dict.fromkeys(atoms.get_chemical_symbols())

    with open(outcar, "r", encoding="utf-8") as f:
        lines = f.readlines()
        for line in lines:
            if "LDAUU" in line and "U (eV)" not in line:
                line_list = line.split()
                for i, symbol in enumerate(map_.keys()):

                    value = float(line_list[-len(map_.keys()) + i])
                    print(value)
                    map_[symbol] = value

    u_vals = []
    for symbol in atoms.get_chemical_symbols():
        u_vals.append(map_[symbol])

    return np.array(u_vals)


def fetch_configs_vasp(calc_dirs: list[str]) -> list[ase.Atoms]:
    """
    from a set of directories containing vasp in files, read configs.
    """
    configs = []
    for directory in calc_dirs:
        if havePOSCAR(directory):
            atoms = ase.io.read(f"{directory}/POSCAR")
            configs.append(atoms)
        else:
            print(f"No POSCAR found in directory: {directory}")

    return configs


# -------------------set ICHARG across database------------------


def set_icharg(value: int, vasp_calc_dir_param: str):

    with open(f"{vasp_calc_dir_param}/INCAR", "r", encoding="utf-8") as f:
        incar_lines = f.readlines()
    found = False
    for i, line in enumerate(incar_lines):
        if "ICHARG" in line:
            incar_lines[i] = f"ICHARG = {value}\n"
            found = True
    if not found:
        incar_lines.append("\n")
        incar_lines.append(f"ICHARG = {value}\n")

    with open(f"{vasp_calc_dir_param}/INCAR", "w", encoding="utf-8") as f:
        new_file_str = "".join(incar_lines)
        f.write(new_file_str)


# -------------------MAGMOM for large databases------------------


def set_magmom(
    magmom_motif_config: ase.Atoms = None,
    vasp_calc_dirs: str = "./QM_calculations",
    atol: float = 1e-4,
    magmom_motif_multi_samples: list[ase.Atoms] = None,
) -> None:
    """
    Given a set of vasp calculation directories, the supercell size, a motif and the magnet moments for the motif, POSCAR is used to set the MAGMOM string.
    Allowing the user to access magnetically ordered states for larger supercells.

    This is a solid specific functionality where

    Parameters
    ----------
    supercell_size :class:`np.ndarray`
        3D array defining supercell size
    motif: :class:`np.ndarray`
        motif of a relaxed solid structure.
    magmom_motif: :class:`np.ndarray`
        magnetic moments of the motif, order of magmom_motif must equal the order of motif. i.e. the magnetic moment of atom located at motif[i] is magmom_motif[i].
    magmom_motif_multi_samples:
        used to define a list of magnetic moment configs used in sampling.

    Returns
    -------
    None : None
        edits INCAR files in call sub directories.
    """

    path = Path(vasp_calc_dirs)
    subdirs = [p for p in path.iterdir() if p.is_dir()]
    for vasp_calc in subdirs:
        if haveINCAR(str(vasp_calc)) and havePOSCAR(str(vasp_calc)):
            if magmom_motif_multi_samples:
                # randomly select a magnetic order to sample
                selection = np.random.choice(len(magmom_motif_multi_samples))
                set_magmom_one_directory(
                    magmom_motif_multi_samples[selection], str(vasp_calc), atol
                )
            else:
                set_magmom_one_directory(magmom_motif_config, str(vasp_calc), atol)

    print(f"Magnetic Moments updated in all vasp sub directories of {vasp_calc_dirs}")


def set_magmom_one_directory(
    magmom_motif_config: ase.Atoms, vasp_calc_dir_param: str, atol: float = 1e-4
) -> None:
    """
    Called on each directory by set_magmom
    """
    # This function is quite brute force and is oppitunity to optimize.

    # read data
    minimal_cell = np.array(magmom_motif_config.cell)
    try:
        motif_magmoms = magmom_motif_config.get_initial_magnetic_moments()
        motif_scaled_pos = magmom_motif_config.get_scaled_positions()
    except Exception as exc:
        raise ValueError(
            "magnetic moment motif input file must contain positions and initial moments."
        ) from exc

    input_atoms = ase.io.read(f"{vasp_calc_dir_param}/POSCAR")
    cell = np.array(input_atoms.cell)

    # create a high symmetry configuration for an expanded cell
    expanded_magmom_motif_config = magmom_motif_config.copy()

    scaling_factor = round(np.linalg.det(cell) / np.linalg.det(minimal_cell))
    volume_scaling_factor = (np.linalg.det(cell) / scaling_factor) / np.linalg.det(
        minimal_cell
    )
    expanded_cell = (volume_scaling_factor) ** (1 / 3) * minimal_cell

    expanded_magmom_motif_config.set_cell(expanded_cell)
    expanded_magmom_motif_config.set_scaled_positions(motif_scaled_pos)
    expanded_magmom_motif_config.set_initial_magnetic_moments(motif_magmoms)

    # now expand this cell to the input supercell.

    S = utils.get_supercell_matrix(cell, expanded_cell)

    if np.allclose(S, np.round(S), atol=atol):
        supercell_magmom_motif = ase.build.make_supercell(
            expanded_magmom_motif_config, np.round(S)
        )
    else:
        raise ValueError(
            "cannot set magmom for configurations with non-interger supercell matricies, maybe need to use a conversion before using as input data."
        )

    final_cell = supercell_magmom_motif.cell
    final_elements = supercell_magmom_motif.get_chemical_symbols()
    final_magmoms = supercell_magmom_motif.get_initial_magnetic_moments()

    possible_vectors = []
    # by expanding range to (-1,1), variations of wrapped co-ords outputed by the MD calculation are dealt with.
    for i, j, k in product(range(-1, 2), range(-1, 2), range(-1, 2)):
        possible_vectors.append(np.array([i, j, k]))
    expected_positions = []  # expected for a relaxed lattice
    expected_mag_moments = []
    expected_elements = []
    for vecs in possible_vectors:
        for l, equilibrium_pos in enumerate(supercell_magmom_motif.get_positions()):
            pos = equilibrium_pos + vecs @ final_cell
            expected_positions.append((pos))
            expected_elements.append(final_elements[l])
            expected_mag_moments.append(final_magmoms[l])  # set corresponding magmom

    # Now match magmom according to minimum distance and element.

    A = input_atoms.positions
    B = np.array(expected_positions)
    all_diff = B[None, :, :] - A[:, None, :]
    all_dist2 = np.sum(all_diff**2, axis=2)

    magmoms = np.zeros((len(input_atoms), 3))

    for i, atom in enumerate(input_atoms):
        symbol = atom.symbol
        element_indices = [
            index for index, val in enumerate(expected_elements) if val == symbol
        ]
        B_this_element = B[np.array(element_indices), :]
        diff = B_this_element - A[i]
        dist2 = np.sum(diff**2, axis=1)
        closest_index = np.argmin(dist2)
        true_index = np.where(dist2[closest_index] == all_dist2)[1][0]
        magmoms[i] = expected_mag_moments[true_index]

    # with magmom known can define magmom str
    magmom_str = "MAGMOM = "
    for i in range(len(input_atoms.positions)):
        mx, my, mz = magmoms[i][0:3]
        magmom_str += f"{mx} {my} {mz} "

    input_atoms.set_initial_magnetic_moments(magmoms)
    writeMAGMOM(f"{vasp_calc_dir_param}/INCAR", new_magmom_str=magmom_str)


def writeMAGMOM(incar: str, new_magmom_str: str) -> None:
    """
    given the path to an INCAR file, writes or updates the MAGMOM string
    """

    with open(incar, "r", encoding="utf-8") as f:
        incar_lines = f.readlines()
    found = False
    for i, line in enumerate(incar_lines):
        if "MAGMOM" in line:
            incar_lines[i] = new_magmom_str + "\n"
            found = True
    if not found:
        incar_lines.append("\n")
        incar_lines.append(new_magmom_str + "\n")

    with open(incar, "w", encoding="utf-8") as f:
        new_file_str = "".join(incar_lines)
        f.write(new_file_str)


# -------------------KPOINTS for large databases------------------


def set_kpoints(
    kspacing: float, grid_type: str = "Gamma", vasp_calc_dirs: str = "./QM_calculations"
):

    path = Path(vasp_calc_dirs)
    subdirs = [p for p in path.iterdir() if p.is_dir()]
    for vasp_calc in subdirs:
        if havePOSCAR(str(vasp_calc)):
            set_kpoints_one_directory(str(vasp_calc), kspacing, grid_type)


def set_kpoints_one_directory(
    vasp_calc_dir_param: str, kspacing: float, grid_type: str = "Gamma"
):
    """

    :param kspacing: units: 2pi/A

    """

    input_atoms = ase.io.read(f"{vasp_calc_dir_param}/POSCAR")
    rcp_lattice_vectors = (
        input_atoms.get_reciprocal_cell() * 2 * np.pi
    )  # reciprocal_lattice_vectors(cell)

    kpoints = np.zeros(3)
    for i in range(3):
        # There are some subtlties
        kpoints[i] = max(
            1, round(np.linalg.norm(rcp_lattice_vectors[i]) / (2 * np.pi * kspacing))
        )

    with open(f"{vasp_calc_dir_param}/KPOINTS", "w", encoding="utf-8") as f:

        f.write(f"K-Spacing Value to Generate K-Mesh: {kspacing:.3f}\n")
        f.write("0\n")
        f.write(f"{grid_type}\n")
        f.write(f"{kpoints[0]:.0f} {kpoints[1]:.0f} {kpoints[2]:.0f}\n")
        f.write("0.0 0.0 0.0\n")


def reciprocal_lattice_vectors(lattice_vectors: np.ndarray):
    """
    Returns the reciprocal lattice vectors of 3x3 lattice vectors.
    """

    rcp_lattice_vectors = np.zeros_like(lattice_vectors)
    V = np.dot(lattice_vectors[0], np.cross(lattice_vectors[1], lattice_vectors[2]))

    pertubations = [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
    for p in pertubations:
        rcp_lattice_vectors[p[0]] = (
            2 * np.pi / V * (np.cross(lattice_vectors[p[1]], lattice_vectors[p[2]]))
        )

    return rcp_lattice_vectors


def set_Uvalue(vasp_calc_dir_param: str, u_vals_param: list[float]):

    if haveINCAR(vasp_calc_dir_param):

        with open(vasp_calc_dir_param + "/INCAR", "r", encoding="utf-8") as f:
            lines = f.readlines()
        updated_lines = []
        for line in lines:
            if "LDAUU" in line:
                new_str = (
                    "LDAUU = " + " ".join([f"{i:.3f}" for i in u_vals_param]) + "\n"
                )
                updated_lines.append(new_str)
            else:
                updated_lines.append(line)

        with open(vasp_calc_dir_param + "/INCAR", "w", encoding="utf-8") as f:
            for line in updated_lines:
                f.write(line)

    else:
        raise ValueError("Tried to write new U value but no INCAR found")


# -------------------ANY INCAR PARAM for large databases------------------

# want a way to set an INCAR parameter given a cell type condition, e.g. number of atoms.


def haveINCAR(directory: str):
    """
    checks if a directory contains INCAR
    """
    path = Path(directory)
    files = [str(p.name) for p in path.iterdir()]
    return "INCAR" in files


def havePOSCAR(directory: str):
    """
    checks if a directory contains POSCAR
    """
    path = Path(directory)
    files = [str(p.name) for p in path.iterdir()]
    return "POSCAR" in files


def haveKPOINTS(directory: str):
    """
    checks if a directory contains KPOINTS
    """
    path = Path(directory)
    files = [str(p.name) for p in path.iterdir()]
    return "KPOINTS" in files


if __name__ == "__main__":

    vasp_calc_dir_arg = sys.argv[1]
    option = sys.argv[2]

    if option == "set_LDAUU":

        num_params = len(sys.argv)
        u_vals_arg = [float(sys.argv[i]) for i in range(3, num_params)]
        set_Uvalue(vasp_calc_dir_arg, u_vals_arg)

    else:
        raise ValueError(f"No option to run vasp.py with {option}")
