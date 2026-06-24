"""Parse and validate multi-system initial structures manifest for active learning."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional
import os

import yaml


@dataclass
class SystemManifestEntry:
    name: str
    reference_pdb: str
    simulation_config: str
    reference_sdf: Optional[str] = None
    reference_xyz: Optional[str] = None
    source_structure: Optional[str] = None


def _abs_path(path: str, base_dir: Optional[str] = None) -> str:
    if base_dir is not None and not os.path.isabs(path):
        return os.path.abspath(os.path.join(base_dir, path))
    return os.path.abspath(path)


def _validate_simulation_yaml(sim_yaml: dict, name: str) -> None:
    if "System" not in sim_yaml:
        raise ValueError(f"simulation_config for {name} must have 'System' section")
    if "Simulation" not in sim_yaml:
        raise ValueError(f"simulation_config for {name} must have 'Simulation' section")
    plumed_config = (
        sim_yaml.get("Simulation", {})
        .get("sampling", {})
        .get("params", {})
        .get("plumed_config")
    )
    if plumed_config is None:
        raise ValueError(
            f"simulation_config for {name} must have "
            "Simulation.sampling.params.plumed_config"
        )


def load_systems_manifest(manifest_path: str) -> List[SystemManifestEntry]:
    manifest_path = os.path.abspath(manifest_path)
    manifest_dir = os.path.dirname(manifest_path)
    with open(manifest_path, "r") as handle:
        data = yaml.load(handle, Loader=yaml.FullLoader)

    systems = data.get("systems")
    if not systems:
        raise ValueError("Manifest must contain a non-empty 'systems' list")

    entries: List[SystemManifestEntry] = []
    for idx, system in enumerate(systems):
        if not isinstance(system, dict):
            raise ValueError(f"systems[{idx}] must be a mapping")

        name = system.get("name")
        if not name:
            raise ValueError(f"systems[{idx}] must have 'name'")

        simulation_config = system.get("simulation_config")
        if not simulation_config:
            raise ValueError(f"systems[{idx}] ({name}) must have 'simulation_config'")
        simulation_config = _abs_path(simulation_config, manifest_dir)
        if not os.path.exists(simulation_config):
            raise FileNotFoundError(
                f"simulation_config not found for {name}: {simulation_config}"
            )

        with open(simulation_config, "r") as handle:
            sim_yaml = yaml.load(handle, Loader=yaml.FullLoader)
        _validate_simulation_yaml(sim_yaml, name)

        reference_pdb = system.get("reference_pdb")
        if not reference_pdb:
            raise ValueError(f"systems[{idx}] ({name}) must have 'reference_pdb'")
        reference_pdb = _abs_path(reference_pdb, manifest_dir)
        if not os.path.exists(reference_pdb):
            raise FileNotFoundError(f"reference_pdb not found for {name}: {reference_pdb}")

        reference_sdf = system.get("reference_sdf")
        if reference_sdf:
            reference_sdf = _abs_path(reference_sdf, manifest_dir)
            if not os.path.exists(reference_sdf):
                raise FileNotFoundError(
                    f"reference_sdf not found for {name}: {reference_sdf}"
                )

        reference_xyz = system.get("reference_xyz")
        if reference_xyz:
            reference_xyz = _abs_path(reference_xyz, manifest_dir)
            if not os.path.exists(reference_xyz):
                raise FileNotFoundError(
                    f"reference_xyz not found for {name}: {reference_xyz}"
                )
        else:
            structure_file = sim_yaml.get("System", {}).get("structure_file")
            if not structure_file:
                raise ValueError(
                    f"systems[{idx}] ({name}): provide 'reference_xyz' or set "
                    "System.structure_file in simulation_config"
                )
            reference_xyz = _abs_path(
                structure_file, os.path.dirname(simulation_config)
            )
            if not os.path.exists(reference_xyz):
                raise FileNotFoundError(
                    f"structure_file not found for {name}: {reference_xyz}"
                )

        plumed_config = (
            sim_yaml.get("Simulation", {})
            .get("sampling", {})
            .get("params", {})
            .get("plumed_config", {})
        )
        yaml_pdb = plumed_config.get("reference_pdb_file")
        if yaml_pdb:
            yaml_pdb_abs = _abs_path(yaml_pdb, os.path.dirname(simulation_config))
            if os.path.normpath(yaml_pdb_abs) != os.path.normpath(reference_pdb):
                raise ValueError(
                    f"systems[{idx}] ({name}): reference_pdb in manifest "
                    f"({reference_pdb}) does not match "
                    f"plumed_config.reference_pdb_file ({yaml_pdb_abs}) "
                    "in simulation_config"
                )

        if plumed_config.get("proton_transfer"):
            raise NotImplementedError(
                "Multi-system active learning does not support proton_transfer "
                f"for system {name}"
            )

        entries.append(
            SystemManifestEntry(
                name=name,
                reference_pdb=reference_pdb,
                simulation_config=simulation_config,
                reference_sdf=reference_sdf,
                reference_xyz=reference_xyz,
                source_structure=reference_xyz,
            )
        )

    return entries
