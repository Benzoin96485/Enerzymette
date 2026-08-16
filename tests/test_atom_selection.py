"""Tests for atom_selection helpers."""

from __future__ import annotations

import pytest

from enerzymette.plumed_config_generator.atom_selection import (
    AtomSpec,
    BondPairSpec,
    atom_spec_from_mapping,
    bond_pair_from_mapping,
    coerce_bond_pairs,
    parse_pdb_atoms,
    resolve_atom_index,
    resolve_bond_pair_indices,
    to_ase_index,
    to_plumed_index,
)


def _write_pdb(path, lines):
    path.write_text("\n".join(lines) + "\n")


def test_atom_spec_rejects_mixed_and_empty():
    with pytest.raises(ValueError, match="cannot mix"):
        AtomSpec(index=1, atom_name="SD", residue_name="SAM").validate()
    with pytest.raises(ValueError, match="either index or PDB"):
        AtomSpec().validate()
    with pytest.raises(ValueError, match="atom_name is required"):
        AtomSpec(residue_name="SAM").validate()
    with pytest.raises(ValueError, match="residue_name and/or residue_number"):
        AtomSpec(atom_name="SD").validate()


def test_atom_spec_from_mapping_aliases(tmp_path):
    spec = atom_spec_from_mapping(
        {"chain": "A", "resname": "SAM", "resid": 10, "atomname": "SD"},
        label="sulphur",
    )
    assert spec == AtomSpec(
        chain_id="A",
        residue_name="SAM",
        residue_number=10,
        atom_name="SD",
    )
    with pytest.raises(ValueError, match="conflicting aliases"):
        atom_spec_from_mapping({"resname": "SAM", "residue_name": "SAM", "atom_name": "SD"})
    with pytest.raises(ValueError, match="unknown fields"):
        atom_spec_from_mapping({"index": 1, "not_a_field": 2})


def test_bond_pair_requires_both_atoms():
    with pytest.raises(ValueError, match="both atom1 and atom2"):
        bond_pair_from_mapping({"atom1": {"index": 1}}, label="forming_bond")


def test_parse_and_resolve_pdb(tmp_path):
    pdb = tmp_path / "toy.pdb"
    _write_pdb(
        pdb,
        [
            # columns roughly matching PDB: atom name 13-16, resname 18-20, chain 22, resseq 23-26
            "HETATM    1  SD  SAM A  10       0.000   0.000   0.000  1.00  0.00           S",
            "HETATM    2  CE  SAM A  10       1.500   0.000   0.000  1.00  0.00           C",
            "HETATM    3  O3  DHB A  20       3.000   0.000   0.000  1.00  0.00           O",
            "HETATM    4  O3  DHB B  20       4.000   0.000   0.000  1.00  0.00           O",
        ],
    )
    records = parse_pdb_atoms(str(pdb))
    assert len(records) == 4
    assert records[0].atom_name == "SD"
    assert records[0].residue_name == "SAM"
    assert records[0].chain_id == "A"
    assert records[0].residue_number == 10

    sd = resolve_atom_index(
        AtomSpec(residue_name="SAM", atom_name="SD"),
        idx_start_from=1,
        reference_pdb_file=str(pdb),
        label="SD",
    )
    assert sd == 1
    assert to_ase_index(sd, 1) == 0
    assert to_plumed_index(sd, 1) == 1

    # Ambiguous without chain
    with pytest.raises(ValueError, match="ambiguous"):
        resolve_atom_index(
            AtomSpec(residue_name="DHB", atom_name="O3"),
            idx_start_from=0,
            reference_pdb_file=str(pdb),
        )

    assert (
        resolve_atom_index(
            AtomSpec(residue_name="DHB", atom_name="O3"),
            idx_start_from=0,
            reference_pdb_file=str(pdb),
            on_ambiguous="last",
        )
        == 3
    )
    assert (
        resolve_atom_index(
            AtomSpec(residue_name="DHB", atom_name="O3"),
            idx_start_from=0,
            reference_pdb_file=str(pdb),
            on_ambiguous="first",
        )
        == 2
    )

    o3 = resolve_atom_index(
        AtomSpec(chain_id="B", residue_name="DHB", atom_name="O3"),
        idx_start_from=0,
        reference_pdb_file=str(pdb),
    )
    assert o3 == 3

    with pytest.raises(ValueError, match="no PDB atom matched"):
        resolve_atom_index(
            AtomSpec(residue_name="SAM", atom_name="XX"),
            idx_start_from=1,
            reference_pdb_file=str(pdb),
        )


def test_explicit_index_bounds_and_bond_pair():
    with pytest.raises(ValueError, match="out of range"):
        resolve_atom_index(
            AtomSpec(index=5),
            idx_start_from=1,
            n_atoms=3,
            label="atom",
        )
    pair = BondPairSpec(AtomSpec(index=1), AtomSpec(index=2))
    assert resolve_bond_pair_indices(pair, idx_start_from=1, n_atoms=3) == (1, 2)


def test_coerce_bond_pairs_accepts_none_one_or_list():
    assert coerce_bond_pairs(None, label="forming_bonds") == []
    assert coerce_bond_pairs([], label="forming_bonds") == []
    single = {"atom1": {"index": 1}, "atom2": {"index": 2}}
    wrapped = coerce_bond_pairs(single, label="forming_bonds")
    assert len(wrapped) == 1
    assert wrapped[0].atom1.index == 1
    listed = coerce_bond_pairs(
        [
            {"atom1": {"index": 1}, "atom2": {"index": 2}},
            BondPairSpec(AtomSpec(index=3), AtomSpec(index=4)),
        ],
        label="forming_bonds",
    )
    assert [pair.atom1.index for pair in listed] == [1, 3]
    with pytest.raises(TypeError, match="expected a bond pair"):
        coerce_bond_pairs("not-a-pair", label="forming_bonds")


def test_pdb_selection_requires_reference_file():
    with pytest.raises(ValueError, match="reference_pdb_file is required"):
        resolve_atom_index(
            AtomSpec(residue_name="SAM", atom_name="SD"),
            idx_start_from=1,
        )
