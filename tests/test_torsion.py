"""Tests for TorsionConfigGenerator."""

from __future__ import annotations

import math

import numpy as np
import pytest
from ase import Atoms

from enerzymette.plumed_config_generator import (
    get_config_generator_class,
    get_config_generator_name,
    list_plumed_cv_plugin_keys,
    resolve_scan_endpoints,
)
from enerzymette.plumed_config_generator.torsion import TorsionConfigGenerator


def _dihedral_atoms(phi_deg: float) -> Atoms:
    """Four atoms with a signed dihedral of approximately ``phi_deg``.

    Bond 1–2 lies on +x. Atom 0 is in the xy plane; atom 3 is rotated about x.
    """
    phi = math.radians(phi_deg)
    return Atoms(
        symbols=["C", "C", "C", "C"],
        positions=[
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, math.cos(phi), math.sin(phi)],
        ],
    )


def test_plugin_registry_includes_torsion():
    assert "torsion" in list_plumed_cv_plugin_keys()
    assert get_config_generator_name("torsion") == "TorsionConfigGenerator"
    assert get_config_generator_class("torsion") is TorsionConfigGenerator


def test_requires_four_atoms():
    system = _dihedral_atoms(0.0)
    with pytest.raises(ValueError, match="atom1, atom2, atom3, and atom4"):
        TorsionConfigGenerator(system, idx_start_from=1, atom1=1, atom2=2)
    with pytest.raises(ValueError, match="Provide atoms"):
        TorsionConfigGenerator(system, idx_start_from=1)


def test_explicit_indices_plumed_numbers_and_scan():
    system = _dihedral_atoms(30.0)
    gen = TorsionConfigGenerator(
        system,
        idx_start_from=1,
        atoms=[1, 2, 3, 4],
    )
    expected = math.radians(system.get_dihedral(0, 1, 2, 3, mic=False))
    assert gen.calc_main_rc() == pytest.approx(expected)
    cv_name, definition = gen.define_main_rc()
    assert cv_name == "phi"
    assert definition == "phi: TORSION ATOMS=1,2,3,4 NOPBC"
    assert "MOLINFO" not in definition
    assert "@" not in definition
    assert gen.get_indices() == {"atom1": 1, "atom2": 2, "atom3": 3, "atom4": 4}

    lines = gen.scan(target_value=0.5, dump_interval=20)
    joined = "\n".join(lines)
    assert "phi: TORSION ATOMS=1,2,3,4 NOPBC" in lines
    assert any(line.startswith("r: RESTRAINT ARG=phi AT=0.5") for line in lines)
    assert "MOLINFO" not in joined
    assert "@phi" not in joined
    assert "@psi" not in joined


def test_idx_start_from_zero():
    system = _dihedral_atoms(0.0)
    gen = TorsionConfigGenerator(
        system,
        idx_start_from=0,
        atom1=0,
        atom2=1,
        atom3=2,
        atom4=3,
    )
    _, definition = gen.define_main_rc()
    assert definition == "phi: TORSION ATOMS=1,2,3,4 NOPBC"


def test_pdb_selectors_emit_numeric_atoms_only(tmp_path):
    pdb = tmp_path / "torsion.pdb"
    pdb.write_text(
        "\n".join(
            [
                "HETATM    1  C1  LIG A   1       0.000   1.000   0.000  1.00  0.00           C",
                "HETATM    2  C2  LIG A   1       0.000   0.000   0.000  1.00  0.00           C",
                "HETATM    3  C3  LIG A   1       1.000   0.000   0.000  1.00  0.00           C",
                "HETATM    4  C4  LIG A   1       1.000   1.000   0.000  1.00  0.00           C",
            ]
        )
        + "\n"
    )
    system = _dihedral_atoms(0.0)
    gen = TorsionConfigGenerator(
        system,
        idx_start_from=1,
        reference_pdb_file=str(pdb),
        atoms=[
            {"resname": "LIG", "atom_name": "C1"},
            {"index": 2},
            {"residue_name": "LIG", "atom_name": "C3"},
            {"chain_id": "A", "residue_name": "LIG", "residue_number": 1, "atom_name": "C4"},
        ],
    )
    _, definition = gen.define_main_rc()
    assert definition == "phi: TORSION ATOMS=1,2,3,4 NOPBC"
    assert "MOLINFO" not in definition
    assert "@" not in definition
    assert gen.get_indices()["atom4"] == 4


def test_resolve_scan_endpoints_full_turn_is_two_pi():
    system = _dihedral_atoms(0.0)
    x0, x1, num, rc = resolve_scan_endpoints(
        system,
        1,
        "torsion",
        {
            "atoms": [1, 2, 3, 4],
            "lower_bound": -math.pi,
            "upper_bound": math.pi,
            "dump_interval": 10,
        },
        num=25,
    )
    assert x0 == pytest.approx(-math.pi)
    assert x1 == pytest.approx(math.pi)
    assert x1 - x0 == pytest.approx(2 * math.pi)
    assert num == 25
    assert rc.cv_name == "phi"
    # Enerzyme interpolates this interval literally.
    grid = np.linspace(x0, x1, num)
    assert grid[-1] - grid[0] == pytest.approx(2 * math.pi)
    assert not np.allclose(grid, grid[0])


def test_resolve_scan_endpoints_target_is_literal_not_shortest_arc():
    system = _dihedral_atoms(170.0)
    current = math.radians(system.get_dihedral(0, 1, 2, 3, mic=False))
    x0, x1, num, rc = resolve_scan_endpoints(
        system,
        1,
        "torsion",
        {
            "atoms": [1, 2, 3, 4],
            "lower_bound": -math.pi,
            "upper_bound": math.pi,
            "dump_interval": 10,
        },
        num=25,
        target_value=-math.pi,
    )
    assert x0 == pytest.approx(current)
    assert x1 == pytest.approx(-math.pi)
    # Literal path ~ −2π, not the ~20° shortest arc.
    assert x1 - x0 == pytest.approx(-math.pi - current)
    assert abs(x1 - x0) > math.pi
    periodic_arc = abs((x1 - x0 + math.pi) % (2 * math.pi) - math.pi)
    assert abs(x1 - x0) > periodic_arc + 1.0
    assert rc.initial_value == pytest.approx(x0)


def test_angle_unit_deg_writes_radians():
    system = _dihedral_atoms(0.0)
    gen = TorsionConfigGenerator(
        system,
        idx_start_from=1,
        atoms=[1, 2, 3, 4],
        angle_unit="deg",
        lower_bound=-180,
        upper_bound=180,
        dump_interval=10,
    )
    rc = gen.build_reaction_coordinate(lower_bound=-180, upper_bound=180, dump_interval=10)
    assert rc.lower_bound == pytest.approx(-math.pi)
    assert rc.upper_bound == pytest.approx(math.pi)

    x0, x1, _, _ = resolve_scan_endpoints(
        system,
        1,
        "torsion",
        {
            "atoms": [1, 2, 3, 4],
            "angle_unit": "deg",
            "lower_bound": -180,
            "upper_bound": 180,
            "dump_interval": 10,
        },
        num=5,
    )
    assert x0 == pytest.approx(-math.pi)
    assert x1 == pytest.approx(math.pi)

    # Enerzyme scan points are radians; AT= must not be left in degrees.
    lines = gen.scan(target_value=math.pi, lower_bound=-180, upper_bound=180, dump_interval=10)
    restraint = next(line for line in lines if line.startswith("r: RESTRAINT ARG=phi"))
    assert "AT=180" not in restraint
    assert "AT=3.141" in restraint or "AT=3.14" in restraint


def test_default_bounds_are_plus_minus_pi():
    system = _dihedral_atoms(45.0)
    gen = TorsionConfigGenerator(system, idx_start_from=1, atoms=[1, 2, 3, 4])
    rc = gen.build_reaction_coordinate()
    assert rc.lower_bound == pytest.approx(-math.pi)
    assert rc.upper_bound == pytest.approx(math.pi)
    x0, x1, _, _ = resolve_scan_endpoints(
        system,
        1,
        "torsion",
        {"atoms": [1, 2, 3, 4], "dump_interval": 5},
        num=3,
    )
    assert x0 == pytest.approx(-math.pi)
    assert x1 == pytest.approx(math.pi)


def test_custom_cv_name():
    system = _dihedral_atoms(0.0)
    gen = TorsionConfigGenerator(
        system,
        idx_start_from=1,
        atoms=[1, 2, 3, 4],
        cv_name="omega",
    )
    cv_name, definition = gen.define_main_rc()
    assert cv_name == "omega"
    assert definition.startswith("omega: TORSION")
    lines = gen.scan(target_value=0.0)
    assert any("RESTRAINT ARG=omega" in line for line in lines)


def test_apply_plumed_scan_yaml_target_value_deg(tmp_path):
    import ase.io
    from enerzymette.scantoolkit.workflow import apply_plumed_scan_sampling

    xyz = tmp_path / "torsion.xyz"
    atoms = _dihedral_atoms(30.0)
    ase.io.write(xyz, atoms)
    config = {
        "Simulation": {
            "idx_start_from": 1,
            "optimize": {"optimizer": "LBFGS"},
        }
    }
    apply_plumed_scan_sampling(
        config,
        structure_path=str(xyz),
        plumed_patch_key="torsion",
        plumed_cv_config={
            "atoms": [1, 2, 3, 4],
            "angle_unit": "deg",
            "lower_bound": -180,
            "upper_bound": 180,
            "target_value": 0,
            "dump_interval": 10,
        },
        n_steps=5,
    )
    params = config["Simulation"]["sampling"]["params"]
    assert params["x0"] == pytest.approx(math.radians(30.0), abs=0.05)
    assert params["x1"] == pytest.approx(0.0)
    assert "target_value" not in params["plumed_config"]


def test_rejects_mixing_atoms_list_and_named():
    system = _dihedral_atoms(0.0)
    with pytest.raises(ValueError, match="not both"):
        TorsionConfigGenerator(
            system,
            idx_start_from=1,
            atoms=[1, 2, 3, 4],
            atom1=1,
        )
