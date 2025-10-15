from typing import Literal, List, Dict, Any, Tuple, Optional
import os
from shutil import copy
from .logger import logger

def write_orca_neb_in(
    neb_in_path: str,
    wrapper_path: str, 
    n_images: int=25,
    restart: bool=False,
    pre_opt: bool=True,
    use_ts: bool=False,
    optimizer: Literal["LBFGS", "BFGS", "VPO", "FIRE"]="LBFGS",
    constraint_freeze_xyz: List[int]=[],
    charge: int=0,
    multiplicity: int=1,
):
    if n_images < 2:
        raise ValueError("n_images must be greater than 2")
    neb_in_prefix = f"""! ExtOpt NEB-CI

%method
ProgExt "{os.path.abspath(wrapper_path)}"
end

%neb
Product "product.xyz"
Opt_Method {optimizer}
{"" if restart else "#"} Restart_ALLXYZFile "neb_MEP.allxyz"
NImages {n_images - 2}
PreOpt {"true" if pre_opt else "false"}
{"" if use_ts else "#"} TS "ts.xyz"
end

"""
    neb_in_constraints = "%geom Constraints\n" + "\n".join([f"{{ C {i-1} C }}" for i in constraint_freeze_xyz]) + "\nend\nend\n"

    neb_in_suffix = f"""
* xyzfile {charge} {multiplicity} reactant.xyz
"""
    neb_in_str = neb_in_prefix + neb_in_constraints + neb_in_suffix
    with open(neb_in_path, "w") as f:
        f.write(neb_in_str)

def read_mep_trj_xyz(mep_path_xyz: str) -> Dict[str, Any]:
    with open(mep_path_xyz, "r") as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())
    n_lines = len(lines)
    n_images = n_lines // (n_atoms + 2)
    xyzblocks = []
    energies = []
    for i in range(n_images):
        xyzblocks.append(lines[i * (n_atoms + 2):(i + 1) * (n_atoms + 2)])
        energies.append(float(xyzblocks[-1][1].split()[-1]))
    return {
        "xyzblocks": xyzblocks,
        "energies": energies,
        "n_images": n_images,
        "n_atoms": n_atoms
    }

def convert_allxyz_to_trj_xyz(allxyz_path: str, trj_xyz_path: str) -> Dict[str, Any]:
    with open(allxyz_path, "r") as f:
        content = f.read()
    xyzblocks = content.split("\n>\n")
    with open(trj_xyz_path, "w") as f:
        f.write("\n".join(xyzblocks))

def get_mep_path_info(elementary_reaction_path: str) -> Dict[str, Any]:
    allxyz_path = os.path.join(elementary_reaction_path, "neb_MEP.allxyz")
    trj_xyz_path = os.path.join(elementary_reaction_path, "neb_MEP_trj.xyz")
    if os.path.exists(trj_xyz_path):
        return read_mep_trj_xyz(trj_xyz_path)
    elif os.path.exists(allxyz_path):
        convert_allxyz_to_trj_xyz(allxyz_path, trj_xyz_path)
        return read_mep_trj_xyz(trj_xyz_path)
    else:
        return None

def read_energy(mep_path_xyz: str) -> Optional[float]:
    if os.path.exists(mep_path_xyz):
        info = read_mep_trj_xyz(mep_path_xyz)
        energy = info["energies"][-1]
    return energy

def make_backup(elementary_reaction_path: str, filename: str):
    copy_index = 0
    while True:
        copy_index += 1
        trial_path = os.path.join(elementary_reaction_path, f"{filename}.{copy_index}")
        if not os.path.exists(trial_path):
            copy(os.path.join(elementary_reaction_path, filename), trial_path)
            logger.info(f"Backup {filename} made to {trial_path}")
            break

def redirect_output(source_stream, dest_stream):
    for line in iter(source_stream.readline, b''):
        dest_stream.write(line.decode("utf-8"))
        dest_stream.flush()

def parse_neb_csv(neb_csv_path: str) -> Dict[str, Any]:
    with open(neb_csv_path, "r") as f:
        csv_lines = f.readlines()
    reactions = []
    if len(csv_lines) > 1:
        for line in csv_lines[1:]:
            line_strip = line.strip()
            if line_strip:
                line_split = line_strip.split(",")
                reactions.append((line_split[0], line_split[1]))
    return reactions
