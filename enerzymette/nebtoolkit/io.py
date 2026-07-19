from typing import Literal, List, Dict, Any, Tuple, Optional
import os
from shutil import copy

import yaml

from ..altoolkit.get_index import get_indices, index_type_map
from ..logger import logger
from ..scantoolkit.io import _abs_path, _charge_from_reference_pdb


def infer_reference_type(reference_path: str) -> str:
    ext = os.path.splitext(reference_path)[1].lower()
    if ext in {".yaml", ".yml"}:
        return "neb_config"
    return "terachem_input"


def parse_neb_config(neb_config_path: str, output_path: str) -> Dict[str, Dict[str, Any]]:
    """Parse YAML neb_config used as an alternative to a TeraChem reference input.

    Expected keys:
      - reference_pdb (required): PDB used to resolve freeze atom indices
      - freeze_index_types (required): keys of index_type_map (e.g. backbone)
      - charge (optional): if set, skip deriving charge from reference PDB
      - reference_sdf (optional): ligand template for charge derivation when charge omitted
      - multiplicity / spinmult (optional, default 1)
    """
    neb_config_path = os.path.abspath(neb_config_path)
    config_dir = os.path.dirname(neb_config_path)
    with open(neb_config_path, "r") as handle:
        data = yaml.load(handle, Loader=yaml.FullLoader)
    if not isinstance(data, dict):
        raise ValueError(f"NEB config must be a YAML mapping: {neb_config_path}")

    reference_pdb = data.get("reference_pdb")
    if not reference_pdb:
        raise ValueError(f"NEB config must have reference_pdb: {neb_config_path}")
    reference_pdb = _abs_path(reference_pdb, config_dir)
    if not os.path.exists(reference_pdb):
        raise FileNotFoundError(f"reference_pdb not found: {reference_pdb}")

    reference_sdf = data.get("reference_sdf")
    if reference_sdf:
        reference_sdf = _abs_path(reference_sdf, config_dir)
        if not os.path.exists(reference_sdf):
            raise FileNotFoundError(f"reference_sdf not found: {reference_sdf}")

    multiplicity = int(data.get("multiplicity", data.get("spinmult", 1)))

    freeze_index_types: List[str] = data.get("freeze_index_types", [])
    if not freeze_index_types:
        raise ValueError(
            f"NEB config must specify freeze_index_types "
            f"(keys of index_type_map: {', '.join(index_type_map)})"
        )
    unknown_types = [t for t in freeze_index_types if t not in index_type_map]
    if unknown_types:
        raise ValueError(
            f"Unknown freeze_index_types: {unknown_types}. "
            f"Valid options: {', '.join(index_type_map)}"
        )
    # 0-based indices (same convention as scan_config); converted for ORCA via idx_start_from
    constraint_freeze_xyz = get_indices(reference_pdb, 0, freeze_index_types)
    logger.info(
        f"Freeze indices from {freeze_index_types} on {reference_pdb}: "
        f"{len(constraint_freeze_xyz)} atoms"
    )

    if data.get("charge") is not None:
        charge = int(data["charge"])
        logger.info(f"Using charge ({charge}) from NEB config {neb_config_path}")
    else:
        topology_dir = os.path.join(os.path.abspath(output_path), "topology")
        charge = _charge_from_reference_pdb(reference_pdb, topology_dir, reference_sdf)

    return {
        "main": {
            "charge": charge,
            "spinmult": multiplicity,
        },
        "constraint_freeze": {
            "xyz": constraint_freeze_xyz,
        },
    }


def write_orca_neb_in(
    neb_in_path: str,
    wrapper_path: str, 
    n_images: int=25,
    restart: bool=False,
    pre_opt: bool=True,
    use_ts: bool=False,
    optimizer: Literal["LBFGS", "BFGS", "VPO", "FIRE"]="LBFGS",
    coordsys: Literal["cartesian", "redundant"]="redundant",
    constraint_freeze_xyz: List[int]=[],
    charge: int=0,
    multiplicity: int=1,
    min_spring_constant: float=0.01,
    max_spring_constant: float=0.1,
    idx_start_from: int=1,
):
    if n_images < 2:
        raise ValueError("n_images must be greater than 2")
    neb_in_prefix = f"""! ExtOpt NEB-CI

%method
ProgExt "{os.path.abspath(wrapper_path)}"
end

%geom
coordsys {coordsys}
end

%neb
Product "product.xyz"
Opt_Method {optimizer}
{"" if restart else "# "}Restart_ALLXYZFile "neb_MEP.allxyz"
NImages {n_images - 2}
PreOpt {"true" if pre_opt else "false"}
{"" if use_ts else "# "}TS "ts.xyz"
SpringConst {min_spring_constant}
SpringConst2 {max_spring_constant}
end

"""
    orca_freeze_indices = [i - idx_start_from for i in constraint_freeze_xyz]
    neb_in_constraints = (
        "%geom Constraints\n"
        + "\n".join([f"{{ C {i} C }}" for i in orca_freeze_indices])
        + "\nend\nend\n"
    )

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
    return copy_index

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
