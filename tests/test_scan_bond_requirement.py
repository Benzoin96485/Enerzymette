"""ASE bond-scan configs must include i0/i1 unless -pp is used."""

from __future__ import annotations

import pytest
from ase import Atoms
from ase.io import write

from enerzymette.scantoolkit.workflow import (
    apply_bond_scan_sampling,
    require_ase_bond_pair,
    write_standalone_scan_config,
)


def test_require_ase_bond_pair_accepts_indices():
    assert require_ase_bond_pair({"bond": {"i0": 2, "i1": 5}}) == (2, 5)


@pytest.mark.parametrize(
    "constraint_scan",
    [None, {}, {"bond": {}}, {"bond": {"i0": 1}}, {"angle": {"i0": 1, "i1": 2}}],
)
def test_require_ase_bond_pair_rejects_empty(constraint_scan):
    with pytest.raises(ValueError, match="pass -pp"):
        require_ase_bond_pair(constraint_scan)


def test_apply_bond_scan_sampling_rejects_empty_bond(tmp_path):
    xyz = tmp_path / "mol.xyz"
    write(xyz, Atoms("CC", positions=[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]))
    config = {"Simulation": {}}
    with pytest.raises(ValueError, match="pass -pp"):
        apply_bond_scan_sampling(
            config,
            structure_path=str(xyz),
            constraint_scan={"bond": {}},
            n_steps=5,
            idx_start_from=1,
        )
    assert "sampling" not in config["Simulation"]


def test_write_standalone_scan_without_bond_requires_plumed(tmp_path):
    xyz = tmp_path / "mol.xyz"
    write(xyz, Atoms("CC", positions=[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]]))
    config_path = tmp_path / "scan.yaml"
    with pytest.raises(ValueError, match="pass -pp"):
        write_standalone_scan_config(
            str(config_path),
            task="scan",
            initial_structure_path=str(xyz),
            charge=0,
            multiplicity=1,
            constraint_freeze_xyz=[],
            idx_start_from=1,
            constraint_scan={"bond": {}},
            n_steps=5,
        )
    assert not config_path.exists()
