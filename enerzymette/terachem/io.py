from pathlib import Path
from typing import List, Optional, TextIO
from shutil import rmtree
import numpy as np
from ase.units import Bohr, Hartree, Debye
from ase.utils import reader, writer

# Made from NWChem interface

@writer
def write_terachem(fd, atoms, params):
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

    fd.write(f"end\n")
    return coordinates_path, scr_path


@writer
def write_xyz(fd, atoms):
    fd.write(f"{len(atoms)}\n")
    fd.write(f"Geometry from ASE\n")
    for atom in atoms:
        fd.write(f"{atom.symbol} {atom.position[0]} {atom.position[1]} {atom.position[2]}\n")


def clean_scr(scr_path):
    if Path(scr_path).exists():
        rmtree(scr_path)


def read_energy(lines: List[str]) -> Optional[float]:
    """Read energy."""
    energy = None
    for line in lines:
        if 'FINAL ENERGY:' in line:
            try:
                energy = float(line.split()[-2])
            except:
                energy = float('nan')
    if energy is not None:
        return energy * Hartree
    return energy


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


def read_grad(lines: List[str]) -> Optional[np.ndarray]:
    """Read gradients."""
    gradients = None
    start_read_grad = False
    read_grad_idx = 0
    for line in lines:
        if 'Gradient units are Hartree/Bohr' in line:
            start_read_grad = True
            gradients = []
            continue
        if read_grad_idx > 2 and "--" in line:
            start_read_grad = False
            read_grad_idx = 0
            continue
        if start_read_grad:
            read_grad_idx += 1
            if read_grad_idx > 2:
                gx, gy, gz = line.split()
                gradients.append([float(gx), float(gy), float(gz)])
    if gradients is not None:
        gradients = -np.array(gradients) * Hartree / Bohr
    return gradients


@reader
def read_terachem_output(fd):
    """ From the TERACHEM output file: Read Energy and dipole moment
    """
    lines = fd.readlines()

    energy = read_energy(lines)
    dipole = read_dipole(lines)
    forces = read_grad(lines)

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
