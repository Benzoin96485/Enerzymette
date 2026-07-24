"""Tests for BondReactionConfigGenerator."""

from __future__ import annotations

import pytest
from ase import Atoms

from enerzymette.plumed_config_generator import (
    get_config_generator_class,
    get_config_generator_name,
    list_plumed_cv_plugin_keys,
    resolve_scan_endpoints,
)
from enerzymette.plumed_config_generator.bond_reaction import BondReactionConfigGenerator


def _three_atom_line(d_break=1.5, d_form=2.5):
    # atoms: 0=break_a, 1=shared, 2=form_b  along x
    return Atoms(
        symbols=["S", "C", "O"],
        positions=[[0.0, 0.0, 0.0], [d_break, 0.0, 0.0], [d_break + d_form, 0.0, 0.0]],
    )


def test_plugin_registry_includes_bond_reaction():
    assert "bond_reaction" in list_plumed_cv_plugin_keys()
    assert get_config_generator_name("bond_reaction") == "BondReactionConfigGenerator"
    assert get_config_generator_class("bond_reaction") is BondReactionConfigGenerator


def test_requires_at_least_one_bond():
    system = _three_atom_line()
    with pytest.raises(ValueError, match="At least one"):
        BondReactionConfigGenerator(system, idx_start_from=1)


def test_forming_only_cv():
    system = _three_atom_line(d_break=1.5, d_form=2.5)
    gen = BondReactionConfigGenerator(
        system,
        idx_start_from=1,
        forming_bond={
            "atom1": {"index": 2},
            "atom2": {"index": 3},
        },
    )
    assert gen.calc_main_rc() == pytest.approx(2.5)
    cv_name, definition = gen.define_main_rc()
    assert cv_name == "rc"
    assert "d1: DISTANCE ATOMS=2,3 NOPBC" in definition
    assert "d0:" not in definition
    assert "COMBINE ARG=d1" in definition
    assert gen.get_indices() == {"forming_a": 2, "forming_b": 3}


def test_breaking_only_cv():
    system = _three_atom_line(d_break=1.5, d_form=2.5)
    gen = BondReactionConfigGenerator(
        system,
        idx_start_from=0,
        breaking_bond={
            "atom1": {"index": 0},
            "atom2": {"index": 1},
        },
    )
    assert gen.calc_main_rc() == pytest.approx(1.5)
    cv_name, definition = gen.define_main_rc()
    assert cv_name == "rc"
    assert "d0: DISTANCE ATOMS=1,2 NOPBC" in definition
    assert "d1:" not in definition


def test_difference_cv_and_scan():
    system = _three_atom_line(d_break=1.5, d_form=2.5)
    gen = BondReactionConfigGenerator(
        system,
        idx_start_from=1,
        forming_bond={"atom1": {"index": 2}, "atom2": {"index": 3}},
        breaking_bond={"atom1": {"index": 1}, "atom2": {"index": 2}},
        max_bond_length=3.0,
        lower_bound=-2.0,
        upper_bound=2.0,
        dump_interval=20,
    )
    assert gen.calc_main_rc() == pytest.approx(1.0)
    cv_name, definition = gen.define_main_rc()
    assert cv_name == "rc"
    assert "rc: COMBINE ARG=d1,d0 COEFFICIENTS=1,-1 PERIODIC=NO" in definition
    assert "dsort: SORT ARG=d0,d1" in definition
    assert "UPPER_WALLS" in definition

    lines = gen.scan(
        target_value=0.5,
        lower_bound=-2.0,
        upper_bound=2.0,
        dump_interval=20,
    )
    assert any(line.startswith("r: RESTRAINT ARG=rc AT=0.5") for line in lines)


def test_pdb_selectors(tmp_path):
    pdb = tmp_path / "pair.pdb"
    pdb.write_text(
        "\n".join(
            [
                "HETATM    1  SD  SAM A  10       0.000   0.000   0.000  1.00  0.00           S",
                "HETATM    2  CE  SAM A  10       1.500   0.000   0.000  1.00  0.00           C",
                "HETATM    3  O3  DHB A  20       4.000   0.000   0.000  1.00  0.00           O",
            ]
        )
        + "\n"
    )
    system = Atoms(
        symbols=["S", "C", "O"],
        positions=[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [4.0, 0.0, 0.0]],
    )
    gen = BondReactionConfigGenerator(
        system,
        idx_start_from=1,
        reference_pdb_file=str(pdb),
        breaking_bond={
            "atom1": {"resname": "SAM", "atom_name": "SD"},
            "atom2": {"resname": "SAM", "atom_name": "CE"},
        },
        forming_bond={
            "atom1": {"resname": "SAM", "atom_name": "CE"},
            "atom2": {
                "chain_id": "A",
                "residue_name": "DHB",
                "residue_number": 20,
                "atom_name": "O3",
            },
        },
    )
    assert gen.calc_main_rc() == pytest.approx(2.5 - 1.5)
    assert gen.get_indices()["breaking_a"] == 1
    assert gen.get_indices()["forming_b"] == 3


def test_resolve_scan_endpoints_bond_reaction():
    system = _three_atom_line(d_break=1.0, d_form=3.0)
    x0, x1, num, rc = resolve_scan_endpoints(
        system,
        1,
        "bond_reaction",
        {
            "forming_bond": {"atom1": {"index": 2}, "atom2": {"index": 3}},
            "breaking_bond": {"atom1": {"index": 1}, "atom2": {"index": 2}},
            "lower_bound": -2.0,
            "upper_bound": 2.0,
            "dump_interval": 10,
        },
        num=5,
    )
    assert x0 == pytest.approx(2.0)
    assert x1 == pytest.approx(-2.0)  # farther bound from +2
    assert num == 5
    assert rc.cv_name == "rc"
