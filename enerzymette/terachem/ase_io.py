from pathlib import Path
from typing import List, Optional, TextIO, Dict, Union
from shutil import rmtree, move
import numpy as np
from ase import Atoms
from ase.units import Bohr, Hartree, Debye
from ase.utils import reader, writer
from .io import read_energy, read_grad

# Made from NWChem interface

@writer
def write_terachem(fd: TextIO, atoms: Atoms, params: Dict[str, Union[str, int, float]], find_mo_flags: Dict[str, bool]):
    # conventional filename: '<name>.inp'
    important_keywords = set()
    for line in params.get("terachemblocks", "").split("\n"):
        clean_line = line.strip()
        if clean_line.startswith("#") or clean_line.startswith("end") or not clean_line:
            continue
        if clean_line.startswith("coordinates "):
            coordinates_path = clean_line.split()[1]
            important_keywords.add("coordinates")
        if clean_line.startswith("charge "):
            important_keywords.add("charge")
        if clean_line.startswith("mult "):
            important_keywords.add("mult")
        if clean_line.startswith("run "):
            important_keywords.add("run")
        if clean_line.startswith("scrdir "):
            scr_path = clean_line.split()[1]
            important_keywords.add("scrdir")
        if clean_line.startswith("guess "):
            important_keywords.add("guess")
        fd.write(f"{clean_line}\n")

    if "charge" not in important_keywords:
        fd.write(f"charge {params.get('charge', 0)}\n")
    if "mult" not in important_keywords:
        fd.write(f"mult {params.get('mult', 1)}\n")
    if "run" not in important_keywords:
        fd.write(f"run {params.get('run', 'gradient')}\n")
    if "coordinates" not in important_keywords:
        coordinates_path = params.get('coordinates', 'geom.xyz')
        fd.write(f"coordinates {coordinates_path}\n")
    if "scrdir" not in important_keywords:
        scr_path = params.get('scrdir', 'scr')
        fd.write(f"scrdir {scr_path}\n")
    if "guess" not in important_keywords:
        if find_mo_flags["c0"]:
            fd.write(f"guess {find_mo_flags['c0']}\n")
        elif find_mo_flags["ca"] and find_mo_flags["cb"]:
            fd.write(f"guess {find_mo_flags['ca']} {find_mo_flags['cb']}\n")

    fd.write(f"end\n")
    return coordinates_path, scr_path


@writer
def write_xyz(fd, atoms):
    fd.write(f"{len(atoms)}\n")
    fd.write(f"Geometry from ASE\n")
    for atom in atoms:
        fd.write(f"{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}\n")


def clean_scr(scr_path: str, keep_mo: bool=True, label: Optional[str]=None):
    find_mo_flags = {
        "c0": '',
        "ca": '',
        "cb": '',
    }
    scr_path_ = Path(scr_path)
    if scr_path_.exists():
        if keep_mo:
            for mo_name in ["c0", "ca", "cb"]:
                mo_path = scr_path_ / mo_name
                if mo_path.exists():
                    new_mo_path = f"{mo_name}_{label}" if label is not None else mo_name
                    find_mo_flags[mo_name] = new_mo_path
                    move(mo_path, scr_path_.parent / new_mo_path)
        rmtree(scr_path)
    return find_mo_flags


def read_dipole(lines: List[str]) -> Optional[np.ndarray]:
    """Read dipole moment.

    Note that the read dipole moment is for the COM frame of reference.
    """
    dipole = None
    for line in lines:
        if 'DIPOLE MOMENT:' in line:
            x, y, z = line.split()[2:5]
            x = float(x[1:-1])
            y = float(y[:-1])
            z = float(z[:-1])
            dipole = np.array([x, y, z])
    if dipole is not None:
        return dipole * Debye  # Return the last match
    return dipole


@reader
def read_terachem_output(fd):
    """ From the TERACHEM output file: Read Energy and dipole moment
    """
    lines = fd.readlines()

    energy = read_energy(lines)
    if energy is not None:
        energy = energy * Hartree

    dipole = read_dipole(lines)
    forces = read_grad(lines)
    if forces is not None:
        forces = -np.array(forces) * Hartree / Bohr

    results = {
        'energy': energy,
        'dipole': dipole,
        'forces': forces,
    }

    return results


def read_terachem_outputs(directory, stdout_path):
    stdout_path = Path(stdout_path)
    results = {}
    results.update(read_terachem_output(stdout_path))
    return results
