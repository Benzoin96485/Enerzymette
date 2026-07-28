"""Shared path and charge helpers for scan/neb YAML config parsing."""

from __future__ import annotations

import os
import subprocess
from typing import Optional

from .logger import logger


def abs_path(path: str, base_dir: str) -> str:
    """Resolve *path* to an absolute path, relative to *base_dir* if needed."""
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


def charge_from_reference_pdb(
    reference_pdb: str,
    topology_dir: str,
    reference_sdf: Optional[str] = None,
) -> int:
    """Derive total formal charge from a reference PDB via ``enerzyme bond``."""
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
