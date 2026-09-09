"""Atom selection helpers for PLUMED reaction-coordinate generators.

An atom may be identified either by an explicit index (in the caller's
``idx_start_from`` convention) or by a PDB selector
``(chain_id, residue_name, residue_number, atom_name)``.  The two styles
must not be mixed on the same atom.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union


@dataclass(frozen=True)
class AtomSpec:
    """Specification of a single atom.

    Exactly one of the following must be provided:

    * ``index``: explicit atom index in the caller's ``idx_start_from``
      convention.
    * PDB fields: at least ``atom_name`` together with enough residue
      identifiers to uniquely match one ATOM/HETATM record.  ``chain_id``,
      ``residue_name``, and ``residue_number`` may be omitted individually
      when uniqueness does not require them, but at least one residue-level
      field (``residue_name`` or ``residue_number``) should normally be set.
    """

    index: Optional[int] = None
    chain_id: Optional[str] = None
    residue_name: Optional[str] = None
    residue_number: Optional[int] = None
    atom_name: Optional[str] = None

    def uses_explicit_index(self) -> bool:
        return self.index is not None

    def uses_pdb_selector(self) -> bool:
        return any(
            value is not None
            for value in (
                self.chain_id,
                self.residue_name,
                self.residue_number,
                self.atom_name,
            )
        )

    def validate(self, *, label: str = "atom") -> None:
        has_index = self.uses_explicit_index()
        has_pdb = self.uses_pdb_selector()
        if has_index and has_pdb:
            raise ValueError(
                f"{label}: cannot mix explicit index with PDB selector fields"
            )
        if not has_index and not has_pdb:
            raise ValueError(
                f"{label}: provide either index or PDB selector fields "
                "(chain_id / residue_name / residue_number / atom_name)"
            )
        if has_pdb and self.atom_name is None:
            raise ValueError(f"{label}: atom_name is required for PDB selection")
        if has_pdb and self.residue_name is None and self.residue_number is None:
            raise ValueError(
                f"{label}: PDB selection requires residue_name and/or residue_number"
            )


@dataclass(frozen=True)
class BondPairSpec:
    """A pair of atoms that define a bond distance."""

    atom1: AtomSpec
    atom2: AtomSpec

    def validate(self, *, label: str = "bond") -> None:
        self.atom1.validate(label=f"{label}.atom1")
        self.atom2.validate(label=f"{label}.atom2")


def atom_spec_from_mapping(
    data: Optional[Mapping[str, Any]],
    *,
    label: str = "atom",
) -> Optional[AtomSpec]:
    """Build an :class:`AtomSpec` from a mapping (e.g. YAML fragment).

    Returns ``None`` when ``data`` is ``None``.  Empty mappings raise.
    Accepted keys (aliases in parentheses):

    * ``index``
    * ``chain_id`` (``chain``, ``chainID``)
    * ``residue_name`` (``resname``, ``residue_name``, ``resn``)
    * ``residue_number`` (``resid``, ``resseq``, ``residue_number``, ``resi``)
    * ``atom_name`` (``name``, ``atom``, ``atomname``)
    """
    if data is None:
        return None
    if not isinstance(data, Mapping):
        raise TypeError(f"{label}: expected a mapping, got {type(data).__name__}")
    if not data:
        raise ValueError(f"{label}: empty mapping is not a valid atom selector")

    remaining: Dict[str, Any] = dict(data)
    index = _pop_alias(remaining, "index", label=label)
    chain_id = _pop_alias(remaining, "chain_id", "chain", "chainID", label=label)
    residue_name = _pop_alias(
        remaining, "residue_name", "resname", "resn", label=label
    )
    residue_number = _pop_alias(
        remaining, "residue_number", "resid", "resseq", "resi", label=label
    )
    atom_name = _pop_alias(
        remaining, "atom_name", "name", "atom", "atomname", label=label
    )
    unknown = sorted(str(key) for key in remaining.keys())
    if unknown:
        raise ValueError(f"{label}: unknown fields {unknown}")

    if index is not None:
        index = int(index)
    if residue_number is not None:
        residue_number = int(residue_number)
    if chain_id is not None:
        chain_id = str(chain_id).strip() or None
    if residue_name is not None:
        residue_name = str(residue_name).strip() or None
    if atom_name is not None:
        atom_name = str(atom_name).strip() or None

    spec = AtomSpec(
        index=index,
        chain_id=chain_id,
        residue_name=residue_name,
        residue_number=residue_number,
        atom_name=atom_name,
    )
    spec.validate(label=label)
    return spec


def bond_pair_from_mapping(
    data: Optional[Mapping[str, Any]],
    *,
    label: str = "bond",
) -> Optional[BondPairSpec]:
    """Build a :class:`BondPairSpec` from a mapping with ``atom1`` / ``atom2``."""
    if data is None:
        return None
    if not isinstance(data, Mapping):
        raise TypeError(f"{label}: expected a mapping, got {type(data).__name__}")
    if not data:
        raise ValueError(f"{label}: empty mapping is not a valid bond pair")

    atom1_raw = data.get("atom1")
    atom2_raw = data.get("atom2")
    if atom1_raw is None or atom2_raw is None:
        raise ValueError(f"{label}: both atom1 and atom2 are required")
    unknown = sorted(str(key) for key in data.keys() if key not in {"atom1", "atom2"})
    if unknown:
        raise ValueError(f"{label}: unknown fields {unknown}")

    atom1 = atom_spec_from_mapping(atom1_raw, label=f"{label}.atom1")
    atom2 = atom_spec_from_mapping(atom2_raw, label=f"{label}.atom2")
    assert atom1 is not None and atom2 is not None
    pair = BondPairSpec(atom1=atom1, atom2=atom2)
    pair.validate(label=label)
    return pair


def _pop_alias(
    data: Dict[str, Any],
    *names: str,
    label: str,
) -> Any:
    """Return the first present alias value; raise if multiple aliases appear."""
    found = [(name, data[name]) for name in names if name in data]
    if len(found) > 1:
        keys = ", ".join(repr(name) for name, _ in found)
        raise ValueError(f"{label}: conflicting aliases {keys}")
    if not found:
        return None
    name, value = found[0]
    data.pop(name)
    return value


@dataclass(frozen=True)
class PdbAtomRecord:
    """One ATOM/HETATM record with a 0-based serial order index."""

    order_index: int
    chain_id: str
    residue_name: str
    residue_number: int
    atom_name: str
    line: str


def parse_pdb_atoms(reference_pdb_file: str) -> Sequence[PdbAtomRecord]:
    """Parse ATOM/HETATM records; ``order_index`` is 0-based encounter order."""
    records = []
    with open(reference_pdb_file, "r") as handle:
        order_index = 0
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            # PDB fixed columns: atom name 12-16, resname 17-20, chain 21,
            # residue number 22-26 (standard PDB; tolerant of common variants).
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain_id = line[21].strip() if len(line) > 21 else ""
            resid_field = line[22:26].strip() if len(line) > 22 else ""
            try:
                residue_number = int(resid_field) if resid_field else 0
            except ValueError as exc:
                raise ValueError(
                    f"Could not parse residue number from PDB line: {line.rstrip()}"
                ) from exc
            records.append(
                PdbAtomRecord(
                    order_index=order_index,
                    chain_id=chain_id,
                    residue_name=residue_name,
                    residue_number=residue_number,
                    atom_name=atom_name,
                    line=line.rstrip("\n"),
                )
            )
            order_index += 1
    return records


def resolve_atom_index(
    spec: AtomSpec,
    *,
    idx_start_from: int,
    reference_pdb_file: Optional[str] = None,
    pdb_records: Optional[Sequence[PdbAtomRecord]] = None,
    n_atoms: Optional[int] = None,
    label: str = "atom",
    on_ambiguous: str = "error",
) -> int:
    """Resolve ``spec`` to an index in the caller's ``idx_start_from`` convention.

    ``on_ambiguous`` controls multi-match PDB selectors:

    * ``"error"`` (default): raise ``ValueError``
    * ``"last"``: keep the last matching atom in file order (historical SAMMT)
    * ``"first"``: keep the first matching atom
    """
    if on_ambiguous not in {"error", "last", "first"}:
        raise ValueError(
            f"on_ambiguous must be 'error', 'last', or 'first'; got {on_ambiguous!r}"
        )
    spec.validate(label=label)

    if spec.uses_explicit_index():
        assert spec.index is not None
        if n_atoms is not None:
            ase_index = spec.index - idx_start_from
            if ase_index < 0 or ase_index >= n_atoms:
                raise ValueError(
                    f"{label}: index {spec.index} is out of range for "
                    f"{n_atoms} atoms with idx_start_from={idx_start_from}"
                )
        return int(spec.index)

    if pdb_records is None:
        if reference_pdb_file is None:
            raise ValueError(
                f"{label}: reference_pdb_file is required for PDB atom selection"
            )
        pdb_records = parse_pdb_atoms(reference_pdb_file)

    matches = [
        record
        for record in pdb_records
        if _pdb_record_matches(record, spec)
    ]
    if not matches:
        raise ValueError(
            f"{label}: no PDB atom matched {_format_pdb_selector(spec)}"
            + (f" in {reference_pdb_file}" if reference_pdb_file else "")
        )
    if len(matches) > 1 and on_ambiguous == "error":
        details = ", ".join(
            f"order={record.order_index} "
            f"{record.chain_id}:{record.residue_name}{record.residue_number}/"
            f"{record.atom_name}"
            for record in matches[:5]
        )
        extra = "" if len(matches) <= 5 else f" (+{len(matches) - 5} more)"
        raise ValueError(
            f"{label}: ambiguous PDB selector {_format_pdb_selector(spec)} "
            f"matched {len(matches)} atoms: {details}{extra}"
        )
    chosen = matches[-1] if on_ambiguous == "last" else matches[0]
    return chosen.order_index + idx_start_from


def resolve_bond_pair_indices(
    pair: BondPairSpec,
    *,
    idx_start_from: int,
    reference_pdb_file: Optional[str] = None,
    pdb_records: Optional[Sequence[PdbAtomRecord]] = None,
    n_atoms: Optional[int] = None,
    label: str = "bond",
) -> Tuple[int, int]:
    """Resolve both atoms of a bond pair to caller-convention indices."""
    pair.validate(label=label)
    if pdb_records is None and reference_pdb_file is not None:
        pdb_records = parse_pdb_atoms(reference_pdb_file)
    index1 = resolve_atom_index(
        pair.atom1,
        idx_start_from=idx_start_from,
        reference_pdb_file=reference_pdb_file,
        pdb_records=pdb_records,
        n_atoms=n_atoms,
        label=f"{label}.atom1",
    )
    index2 = resolve_atom_index(
        pair.atom2,
        idx_start_from=idx_start_from,
        reference_pdb_file=reference_pdb_file,
        pdb_records=pdb_records,
        n_atoms=n_atoms,
        label=f"{label}.atom2",
    )
    return index1, index2


def to_ase_index(index: int, idx_start_from: int) -> int:
    """Convert a caller-convention index to ASE 0-based."""
    return index - idx_start_from


def to_plumed_index(index: int, idx_start_from: int) -> int:
    """Convert a caller-convention index to PLUMED 1-based."""
    return index + 1 - idx_start_from


def _pdb_record_matches(record: PdbAtomRecord, spec: AtomSpec) -> bool:
    if spec.atom_name is not None and record.atom_name != spec.atom_name:
        return False
    if spec.residue_name is not None and record.residue_name != spec.residue_name:
        return False
    if spec.residue_number is not None and record.residue_number != spec.residue_number:
        return False
    if spec.chain_id is not None and record.chain_id != spec.chain_id:
        return False
    return True


def _format_pdb_selector(spec: AtomSpec) -> str:
    parts = []
    if spec.chain_id is not None:
        parts.append(f"chain_id={spec.chain_id!r}")
    if spec.residue_name is not None:
        parts.append(f"residue_name={spec.residue_name!r}")
    if spec.residue_number is not None:
        parts.append(f"residue_number={spec.residue_number}")
    if spec.atom_name is not None:
        parts.append(f"atom_name={spec.atom_name!r}")
    return "{" + ", ".join(parts) + "}"


def coerce_bond_pair(
    value: Optional[Union[BondPairSpec, Mapping[str, Any]]],
    *,
    label: str,
) -> Optional[BondPairSpec]:
    """Accept either a :class:`BondPairSpec` or a YAML-style mapping."""
    if value is None:
        return None
    if isinstance(value, BondPairSpec):
        value.validate(label=label)
        return value
    return bond_pair_from_mapping(value, label=label)


BondPairLike = Union[BondPairSpec, Mapping[str, Any]]


def coerce_bond_pairs(
    value: Optional[Union[BondPairLike, Sequence[BondPairLike]]],
    *,
    label: str,
) -> List[BondPairSpec]:
    """Accept one pair, a sequence of pairs, or ``None``.

    ``None`` and an empty sequence both yield ``[]``.  A single
    :class:`BondPairSpec` or YAML mapping is wrapped as a one-element list.
    """
    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        raise TypeError(
            f"{label}: expected a bond pair mapping, BondPairSpec, or sequence; "
            f"got {type(value).__name__}"
        )
    if isinstance(value, BondPairSpec):
        value.validate(label=label)
        return [value]
    if isinstance(value, Mapping):
        pair = bond_pair_from_mapping(value, label=label)
        if pair is None:
            return []
        return [pair]
    if isinstance(value, (list, tuple)):
        pairs: List[BondPairSpec] = []
        for index, item in enumerate(value):
            pair = coerce_bond_pair(item, label=f"{label}[{index}]")
            if pair is None:
                raise ValueError(f"{label}[{index}]: bond pair cannot be None")
            pairs.append(pair)
        return pairs
    raise TypeError(
        f"{label}: expected a bond pair mapping, BondPairSpec, or sequence; "
        f"got {type(value).__name__}"
    )
