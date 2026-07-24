"""Regression tests for SAMMT as a specialized bond-reaction generator."""

from __future__ import annotations

import pytest
from ase import Atoms

from enerzymette.plumed_config_generator.sammt import (
    SAMMTConfigGenerator,
    get_sammt_index,
    get_sammt_scan_bond_indices,
)
from enerzymette.plumed_config_generator.bond_reaction import BondReactionConfigGenerator


def _sammt_pdb(path):
    path.write_text(
        "\n".join(
            [
                "HETATM    1  SD  SAM A   1       0.000   0.000   0.000  1.00  0.00           S",
                "HETATM    2  CE  SAM A   1       1.800   0.000   0.000  1.00  0.00           C",
                "HETATM    3  O3  DHB A   2       3.500   0.000   0.000  1.00  0.00           O",
            ]
        )
        + "\n"
    )


def test_get_sammt_index_and_scan_bond(tmp_path):
    pdb = tmp_path / "sammt.pdb"
    _sammt_pdb(pdb)
    s, c, n = get_sammt_index(1, str(pdb), "DHB", "O3")
    assert (s, c, n) == (1, 2, 3)
    assert get_sammt_scan_bond_indices(str(pdb), "DHB", "O3") == (1, 2)


def test_get_sammt_index_keeps_last_ambiguous_match(tmp_path):
    """Historical SAMMT scanners overwrote matches; keep the last occurrence."""
    pdb = tmp_path / "dup_sammt.pdb"
    pdb.write_text(
        "\n".join(
            [
                "HETATM    1  SD  SAM A   1       0.000   0.000   0.000  1.00  0.00           S",
                "HETATM    2  CE  SAM A   1       1.800   0.000   0.000  1.00  0.00           C",
                "HETATM    3  O3  DHB A   2       3.500   0.000   0.000  1.00  0.00           O",
                "HETATM    4  SD  SAM B   3       0.000   1.000   0.000  1.00  0.00           S",
                "HETATM    5  CE  SAM B   3       1.800   1.000   0.000  1.00  0.00           C",
                "HETATM    6  O3  DHB B   4       3.500   1.000   0.000  1.00  0.00           O",
            ]
        )
        + "\n"
    )
    assert get_sammt_index(1, str(pdb), "DHB", "O3") == (4, 5, 6)
    assert get_sammt_scan_bond_indices(str(pdb), "DHB", "O3") == (4, 5)

    system = Atoms(
        symbols=["S", "C", "O", "S", "C", "O"],
        positions=[
            [0.0, 0.0, 0.0],
            [1.8, 0.0, 0.0],
            [3.5, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.8, 1.0, 0.0],
            [3.5, 1.0, 0.0],
        ],
    )
    gen = SAMMTConfigGenerator(
        system,
        idx_start_from=1,
        reference_pdb_file=str(pdb),
        substrate="DHB",
        nucleophile="O3",
    )
    assert gen.get_indices()["sulphur"] == 4
    assert gen.get_indices()["methyl_carbon"] == 5
    assert gen.get_indices()["nucleophile"] == 6


def test_sammt_from_pdb_matches_generic_difference(tmp_path):
    pdb = tmp_path / "sammt.pdb"
    _sammt_pdb(pdb)
    system = Atoms(
        symbols=["S", "C", "O"],
        positions=[[0.0, 0.0, 0.0], [1.8, 0.0, 0.0], [3.5, 0.0, 0.0]],
    )
    sammt = SAMMTConfigGenerator(
        system,
        idx_start_from=1,
        reference_pdb_file=str(pdb),
        substrate="DHB",
        nucleophile="O3",
    )
    generic = BondReactionConfigGenerator(
        system,
        idx_start_from=1,
        forming_bond={"atom1": {"index": 2}, "atom2": {"index": 3}},
        breaking_bond={"atom1": {"index": 1}, "atom2": {"index": 2}},
        max_bond_length=3.0,
        cv_name="dd",
    )
    # Force SAMMT-compatible print args on generic for text comparison of CV defs.
    generic.default_print_args = sammt.default_print_args

    assert isinstance(sammt, BondReactionConfigGenerator)
    assert sammt.calc_main_rc() == pytest.approx(generic.calc_main_rc())
    assert sammt.calc_main_rc() == pytest.approx((3.5 - 1.8) - 1.8)

    cv_name, definition = sammt.define_main_rc()
    assert cv_name == "dd"
    assert "d0: DISTANCE ATOMS=1,2 NOPBC" in definition
    assert "d1: DISTANCE ATOMS=2,3 NOPBC" in definition
    assert "dd: COMBINE ARG=d1,d0 COEFFICIENTS=1,-1 PERIODIC=NO" in definition
    assert "UPPER_WALLS" in definition

    assert sammt.get_indices()["sulphur"] == 1
    assert sammt.get_indices()["nucleophile"] == 3
    assert sammt.get_indices()["forming_a"] == 2


def test_sammt_explicit_indices():
    system = Atoms(
        symbols=["S", "C", "O"],
        positions=[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]],
    )
    gen = SAMMTConfigGenerator(
        system,
        idx_start_from=0,
        index_sulphur=0,
        index_methyl_carbon=1,
        index_nucleophile=2,
        max_bond_length=None,
    )
    assert gen.calc_main_rc() == pytest.approx(1.5 - 1.5)
    _, definition = gen.define_main_rc()
    assert "UPPER_WALLS" not in definition
    assert "d0: DISTANCE ATOMS=1,2 NOPBC" in definition


def test_sammt_rejects_generic_bond_kwargs():
    system = Atoms(symbols=["S", "C"], positions=[[0, 0, 0], [1, 0, 0]])
    with pytest.raises(ValueError, match="does not accept forming_bond"):
        SAMMTConfigGenerator(
            system,
            index_sulphur=1,
            index_methyl_carbon=2,
            index_nucleophile=2,
            forming_bond={"atom1": {"index": 1}, "atom2": {"index": 2}},
        )


def test_sammt_steered_md_smoke():
    system = Atoms(
        symbols=["S", "C", "O"],
        positions=[[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [3.0, 0.0, 0.0]],
    )
    gen = SAMMTConfigGenerator(
        system,
        idx_start_from=1,
        index_sulphur=1,
        index_methyl_carbon=2,
        index_nucleophile=3,
        integrate_config={"n_step": 100},
    )
    lines = gen.standard_steered_md(
        lower_bound=-2.0,
        upper_bound=2.0,
        dump_interval=10,
    )
    assert any("MOVINGRESTRAINT ARG=dd" in line for line in lines)
    assert any(line.startswith("PRINT ARG=d0,d1,dsort.1,dd,mr.*") for line in lines)
