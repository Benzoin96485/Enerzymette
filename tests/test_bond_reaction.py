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
    assert gen.get_indices() == {
        "forming_a": 2,
        "forming_b": 3,
        "forming_0_a": 2,
        "forming_0_b": 3,
    }


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
    assert "UPPER_WALLS" not in definition
    assert any("UPPER_WALLS ARG=dsort.1" in line for line in gen.define_additional_restraints())

    lines = gen.scan(
        target_value=0.5,
        lower_bound=-2.0,
        upper_bound=2.0,
        dump_interval=20,
    )
    assert any(line.startswith("r: RESTRAINT ARG=rc AT=0.5") for line in lines)
    assert any("UPPER_WALLS ARG=dsort.1" in line for line in lines)


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


def test_custom_cv_name_in_print_args():
    system = _three_atom_line(d_break=1.5, d_form=2.5)
    gen = BondReactionConfigGenerator(
        system,
        idx_start_from=1,
        forming_bond={"atom1": {"index": 2}, "atom2": {"index": 3}},
        breaking_bond={"atom1": {"index": 1}, "atom2": {"index": 2}},
        cv_name="dd",
        integrate_config={"n_step": 100},
    )
    assert gen.default_print_args == "d0,d1,dsort.1,dd,mr.*"
    lines = gen.standard_steered_md(
        lower_bound=-2.0,
        upper_bound=2.0,
        dump_interval=10,
    )
    assert any("MOVINGRESTRAINT ARG=dd" in line for line in lines)
    assert any(line.startswith("PRINT ARG=d0,d1,dsort.1,dd,mr.*") for line in lines)
    assert not any(",rc," in line or line.endswith(",rc") for line in lines if line.startswith("PRINT"))


def test_forming_cv_name_d1_avoids_undefined_rc():
    system = _three_atom_line()
    gen = BondReactionConfigGenerator(
        system,
        idx_start_from=1,
        forming_bond={"atom1": {"index": 2}, "atom2": {"index": 3}},
        cv_name="d1",
        integrate_config={"n_step": 50},
    )
    assert gen.default_print_args == "d1,mr.*"
    lines = gen.standard_steered_md(
        lower_bound=1.0,
        upper_bound=3.0,
        dump_interval=5,
    )
    print_line = next(line for line in lines if line.startswith("PRINT ARG="))
    assert print_line == "PRINT ARG=d1,mr.* STRIDE=5"
    assert "rc" not in print_line


def _four_atom_line(d_break=1.0, d_form0=2.0, d_form1=2.5):
    return Atoms(
        symbols=["S", "C", "O", "N"],
        positions=[
            [0.0, 0.0, 0.0],
            [d_break, 0.0, 0.0],
            [d_break + d_form0, 0.0, 0.0],
            [d_break + d_form0 + d_form1, 0.0, 0.0],
        ],
    )


def test_rejects_singular_and_plural_together():
    system = _four_atom_line()
    pair = {"atom1": {"index": 2}, "atom2": {"index": 3}}
    with pytest.raises(ValueError, match="not both"):
        BondReactionConfigGenerator(
            system,
            idx_start_from=1,
            forming_bond=pair,
            forming_bonds=[pair],
        )
    with pytest.raises(ValueError, match="not both"):
        BondReactionConfigGenerator(
            system,
            idx_start_from=1,
            breaking_bond={"atom1": {"index": 1}, "atom2": {"index": 2}},
            breaking_bonds=[{"atom1": {"index": 1}, "atom2": {"index": 2}}],
        )


def test_concerted_two_forming_one_breaking():
    system = _four_atom_line(d_break=1.0, d_form0=2.0, d_form1=2.5)
    gen = BondReactionConfigGenerator(
        system,
        idx_start_from=1,
        forming_bonds=[
            {"atom1": {"index": 2}, "atom2": {"index": 3}},
            {"atom1": {"index": 3}, "atom2": {"index": 4}},
        ],
        breaking_bonds=[{"atom1": {"index": 1}, "atom2": {"index": 2}}],
        max_bond_length=4.0,
        lower_bound=-4.0,
        upper_bound=4.0,
        dump_interval=20,
    )
    assert gen.calc_main_rc() == pytest.approx(2.0 + 2.5 - 1.0)
    cv_name, definition = gen.define_main_rc()
    assert cv_name == "rc"
    assert "d0: DISTANCE ATOMS=1,2 NOPBC" in definition
    assert "d1: DISTANCE ATOMS=2,3 NOPBC" in definition
    assert "df1: DISTANCE ATOMS=3,4 NOPBC" in definition
    assert "dsort: SORT ARG=d0,d1,df1" in definition
    assert "rc: COMBINE ARG=d1,df1,d0 COEFFICIENTS=1,1,-1 PERIODIC=NO" in definition
    assert any("UPPER_WALLS ARG=dsort.1" in line for line in gen.define_additional_restraints())
    indices = gen.get_indices()
    assert indices["breaking_a"] == 1
    assert indices["forming_a"] == 2
    assert indices["forming_1_a"] == 3
    assert indices["forming_1_b"] == 4

    lines = gen.scan(
        target_value=0.0,
        lower_bound=-4.0,
        upper_bound=4.0,
        dump_interval=20,
    )
    assert any(line.startswith("r: RESTRAINT ARG=rc AT=0.0") for line in lines)


def test_concerted_two_forming_only():
    system = _four_atom_line(d_break=1.0, d_form0=2.0, d_form1=2.5)
    gen = BondReactionConfigGenerator(
        system,
        idx_start_from=1,
        forming_bonds=[
            {"atom1": {"index": 2}, "atom2": {"index": 3}},
            {"atom1": {"index": 3}, "atom2": {"index": 4}},
        ],
        max_bond_length=4.0,
    )
    assert gen.calc_main_rc() == pytest.approx(2.0 + 2.5)
    cv_name, definition = gen.define_main_rc()
    assert cv_name == "rc"
    assert "d1: DISTANCE ATOMS=2,3 NOPBC" in definition
    assert "df1: DISTANCE ATOMS=3,4 NOPBC" in definition
    assert "d0:" not in definition
    assert "dsort: SORT ARG=d1,df1" in definition
    assert "rc: COMBINE ARG=d1,df1 COEFFICIENTS=1,1 PERIODIC=NO" in definition
    assert any("UPPER_WALLS ARG=dsort.2" in line for line in gen.define_additional_restraints())


def test_resolve_scan_endpoints_concerted_lists():
    system = _four_atom_line(d_break=1.0, d_form0=2.0, d_form1=3.0)
    x0, x1, num, rc = resolve_scan_endpoints(
        system,
        1,
        "bond_reaction",
        {
            "forming_bonds": [
                {"atom1": {"index": 2}, "atom2": {"index": 3}},
                {"atom1": {"index": 3}, "atom2": {"index": 4}},
            ],
            "breaking_bonds": [{"atom1": {"index": 1}, "atom2": {"index": 2}}],
            "lower_bound": -4.0,
            "upper_bound": 4.0,
            "dump_interval": 10,
        },
        num=7,
    )
    assert x0 == pytest.approx(2.0 + 3.0 - 1.0)
    assert x1 == pytest.approx(-4.0)
    assert num == 7
    assert rc.cv_name == "rc"


def test_concerted_pdb_selectors(tmp_path):
    pdb = tmp_path / "concerted.pdb"
    pdb.write_text(
        "\n".join(
            [
                "HETATM    1  SD  SAM A  10       0.000   0.000   0.000  1.00  0.00           S",
                "HETATM    2  CE  SAM A  10       1.000   0.000   0.000  1.00  0.00           C",
                "HETATM    3  SG  CYT A 271       3.000   0.000   0.000  1.00  0.00           S",
                "HETATM    4  O3  DHB A  20       5.500   0.000   0.000  1.00  0.00           O",
            ]
        )
        + "\n"
    )
    system = Atoms(
        symbols=["S", "C", "S", "O"],
        positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [5.5, 0.0, 0.0]],
    )
    gen = BondReactionConfigGenerator(
        system,
        idx_start_from=1,
        reference_pdb_file=str(pdb),
        breaking_bonds=[
            {
                "atom1": {"resname": "SAM", "atom_name": "SD"},
                "atom2": {"resname": "SAM", "atom_name": "CE"},
            }
        ],
        forming_bonds=[
            {
                "atom1": {"residue_name": "SAM", "atom_name": "CE"},
                "atom2": {
                    "chain_id": "A",
                    "residue_name": "CYT",
                    "residue_number": 271,
                    "atom_name": "SG",
                },
            },
            {
                "atom1": {"index": 3},
                "atom2": {"index": 4},
            },
        ],
    )
    assert gen.calc_main_rc() == pytest.approx(2.0 + 2.5 - 1.0)
    indices = gen.get_indices()
    assert indices["breaking_0_a"] == 1
    assert indices["forming_0_b"] == 3
    assert indices["forming_1_b"] == 4
    _, definition = gen.define_main_rc()
    assert "COMBINE ARG=d1,df1,d0 COEFFICIENTS=1,1,-1" in definition
