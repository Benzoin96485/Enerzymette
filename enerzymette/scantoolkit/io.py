import os
import subprocess
from typing import Any, Dict, List, Optional, Tuple

import yaml

from ..altoolkit.get_index import get_indices, index_type_map
from ..logger import logger
from ..plumed_config_generator.sammt import get_sammt_scan_bond_indices
from ..terachem.io import parse_terachem_input, write_terachem_input
import ase.io

_SCAN_BOND_PLUGINS = {
    "sammt": get_sammt_scan_bond_indices,
}


def infer_reference_type(reference_path: str) -> str:
    ext = os.path.splitext(reference_path)[1].lower()
    if ext in {".yaml", ".yml"}:
        return "scan_config"
    return "terachem_input"


def _abs_path(path: str, base_dir: str) -> str:
    if not os.path.isabs(path):
        return os.path.abspath(os.path.join(base_dir, path))
    return os.path.abspath(path)


def _run_enerzyme_bond(
    reference_pdb: str,
    output_img_path: str,
    output_mol_path: str,
    reference_sdf: Optional[str] = None,
) -> None:
    enerzyme_bond_args = [
        "enerzyme", "bond",
        "-p", reference_pdb,
        "-i", output_img_path,
        "-m", output_mol_path,
    ]
    if reference_sdf is not None:
        enerzyme_bond_args.extend(["-t", reference_sdf])
    enerzyme_subprocess = subprocess.Popen(enerzyme_bond_args)
    enerzyme_subprocess.wait()


def _charge_from_reference_pdb(
    reference_pdb: str,
    topology_dir: str,
    reference_sdf: Optional[str] = None,
) -> int:
    from rdkit.Chem import GetFormalCharge, MolFromMolFile

    os.makedirs(topology_dir, exist_ok=True)
    output_img_path = os.path.join(topology_dir, "cluster.png")
    output_mol_path = os.path.join(topology_dir, "cluster.mol")
    _run_enerzyme_bond(
        reference_pdb,
        output_img_path,
        output_mol_path,
        reference_sdf=reference_sdf,
    )
    cluster_mol = MolFromMolFile(output_mol_path, removeHs=False)
    if cluster_mol is None:
        raise ValueError(f"Could not read topology mol written by enerzyme bond: {output_mol_path}")
    charge = GetFormalCharge(cluster_mol)
    logger.info(
        f"Using total charge ({charge}) from reference PDB {reference_pdb}"
        + (f" with template {reference_sdf}" if reference_sdf else "")
    )
    return charge


def _resolve_scan_bond_indices(
    bond_config: dict,
    reference_pdb: str,
) -> Tuple[int, int]:
    if bond_config.get("i0") is not None and bond_config.get("i1") is not None:
        return int(bond_config["i0"]), int(bond_config["i1"])

    plugin = bond_config.get("plugin")
    if plugin is None:
        raise ValueError(
            "constraint_scan.bond must specify plugin (e.g. sammt) or explicit i0/i1"
        )
    resolver = _SCAN_BOND_PLUGINS.get(plugin)
    if resolver is None:
        known = ", ".join(sorted(_SCAN_BOND_PLUGINS))
        raise ValueError(f"Unknown scan bond plugin {plugin!r}. Known plugins: {known}")

    substrate = bond_config.get("substrate")
    nucleophile = bond_config.get("nucleophile")
    if not substrate or not nucleophile:
        raise ValueError(
            f"constraint_scan.bond with plugin={plugin!r} requires substrate and nucleophile"
        )
    pdb_for_indices = bond_config.get("reference_pdb_file", reference_pdb)
    return resolver(pdb_for_indices, substrate, nucleophile)


def parse_scan_config(scan_config_path: str, output_path: str) -> Dict[str, Dict[str, Any]]:
    scan_config_path = os.path.abspath(scan_config_path)
    config_dir = os.path.dirname(scan_config_path)
    with open(scan_config_path, "r") as handle:
        data = yaml.load(handle, Loader=yaml.FullLoader)
    if not isinstance(data, dict):
        raise ValueError(f"Scan config must be a YAML mapping: {scan_config_path}")

    reference_pdb = data.get("reference_pdb")
    if not reference_pdb:
        raise ValueError(f"Scan config must have reference_pdb: {scan_config_path}")
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
            f"Scan config must specify freeze_index_types "
            f"(keys of index_type_map: {', '.join(index_type_map)})"
        )
    unknown_types = [t for t in freeze_index_types if t not in index_type_map]
    if unknown_types:
        raise ValueError(
            f"Unknown freeze_index_types: {unknown_types}. "
            f"Valid options: {', '.join(index_type_map)}"
        )
    constraint_freeze_xyz = get_indices(reference_pdb, 0, freeze_index_types)
    logger.info(
        f"Freeze indices from {freeze_index_types} on {reference_pdb}: "
        f"{len(constraint_freeze_xyz)} atoms"
    )

    constraint_scan_section = data.get("constraint_scan")
    if not isinstance(constraint_scan_section, dict):
        raise ValueError("Scan config must have constraint_scan mapping")
    bond_config = constraint_scan_section.get("bond")
    if not isinstance(bond_config, dict):
        raise ValueError("Scan config constraint_scan must include bond section")

    i0, i1 = _resolve_scan_bond_indices(bond_config, reference_pdb)
    logger.info(f"Scan bond indices from config: i0={i0}, i1={i1}")

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
        "constraint_scan": {
            "bond": {
                "i0": i0,
                "i1": i1,
            },
        },
    }


def update_terachem_scan_input(terachem_input_file: str, updated_structure_xyz: str, output_path: str) -> None:
    info = parse_terachem_input(terachem_input_file)
    i0 = info["constraint_scan"]["bond"]["i0"]
    i1 = info["constraint_scan"]["bond"]["i1"]
    updated_structure = ase.io.read(updated_structure_xyz, index=-1)
    x0 = updated_structure.get_distance(i0 - 1, i1 - 1)
    info["constraint_scan"]["bond"]["x0"] = x0
    write_terachem_input(info, output_path)
