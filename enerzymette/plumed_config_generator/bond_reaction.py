"""Generic forming / breaking bond reaction-coordinate generator.

Users supply one or more complete atom pairs:

* ``forming_bond`` / ``forming_bonds``: bond(s) being formed
* ``breaking_bond`` / ``breaking_bonds``: bond(s) being broken

The main collective variable is chosen automatically:

* forming only  → ``rc = sum(d_forming)``
* breaking only → ``rc = sum(d_breaking)``
* both         → ``rc = sum(d_forming) - sum(d_breaking)``

Singular ``forming_bond`` / ``breaking_bond`` are aliases for a one-element
list.  Do not mix singular and plural keys for the same role.

Atoms may be given as explicit indices or as PDB selectors; see
:mod:`enerzymette.plumed_config_generator.atom_selection`.
"""

from __future__ import annotations

from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

from ase import Atoms
from ase.units import kcal, mol

from enerzymette.plumed_config_generator._engine import PlumedConfigGenerator
from enerzymette.plumed_config_generator.atom_selection import (
    BondPairSpec,
    coerce_bond_pairs,
    resolve_bond_pair_indices,
    to_ase_index,
    to_plumed_index,
)

BondPairLike = Union[BondPairSpec, Mapping[str, Any]]
BondPairInput = Optional[BondPairLike]
BondPairsInput = Optional[Union[BondPairLike, Sequence[BondPairLike]]]


class BondReactionConfigGenerator(PlumedConfigGenerator):
    """PLUMED generator for forming / breaking / difference bond CVs."""

    default_cv_name: str = "rc"
    default_print_args: Optional[str] = None

    def __init__(
        self,
        system: Atoms,
        *,
        forming_bond: BondPairInput = None,
        breaking_bond: BondPairInput = None,
        forming_bonds: BondPairsInput = None,
        breaking_bonds: BondPairsInput = None,
        max_bond_length: Optional[float] = None,
        cv_name: Optional[str] = None,
        **kwargs,
    ) -> None:
        super().__init__(system, **kwargs)
        if forming_bond is not None and forming_bonds is not None:
            raise ValueError("Provide forming_bond or forming_bonds, not both")
        if breaking_bond is not None and breaking_bonds is not None:
            raise ValueError("Provide breaking_bond or breaking_bonds, not both")

        forming_source = forming_bonds if forming_bonds is not None else forming_bond
        breaking_source = breaking_bonds if breaking_bonds is not None else breaking_bond
        self.forming_bonds = coerce_bond_pairs(forming_source, label="forming_bonds")
        self.breaking_bonds = coerce_bond_pairs(breaking_source, label="breaking_bonds")
        if not self.forming_bonds and not self.breaking_bonds:
            raise ValueError(
                "At least one of forming_bond, forming_bonds, breaking_bond, or "
                "breaking_bonds must be provided"
            )
        self.forming_bond = self.forming_bonds[0] if self.forming_bonds else None
        self.breaking_bond = self.breaking_bonds[0] if self.breaking_bonds else None
        self.max_bond_length = max_bond_length
        if cv_name is not None:
            self.default_cv_name = cv_name

        n_atoms = len(system)
        pdb_file = self.reference_pdb
        self.forming_indices: List[Tuple[int, int]] = [
            resolve_bond_pair_indices(
                pair,
                idx_start_from=self.idx_start_from,
                reference_pdb_file=pdb_file,
                n_atoms=n_atoms,
                label=f"forming_bonds[{i}]",
            )
            for i, pair in enumerate(self.forming_bonds)
        ]
        self.breaking_indices: List[Tuple[int, int]] = [
            resolve_bond_pair_indices(
                pair,
                idx_start_from=self.idx_start_from,
                reference_pdb_file=pdb_file,
                n_atoms=n_atoms,
                label=f"breaking_bonds[{i}]",
            )
            for i, pair in enumerate(self.breaking_bonds)
        ]

        self.index_forming_a: Optional[int] = None
        self.index_forming_b: Optional[int] = None
        self.index_breaking_a: Optional[int] = None
        self.index_breaking_b: Optional[int] = None
        if self.forming_indices:
            self.index_forming_a, self.index_forming_b = self.forming_indices[0]
        if self.breaking_indices:
            self.index_breaking_a, self.index_breaking_b = self.breaking_indices[0]

        self._forming_labels = self._distance_labels_for(self.forming_indices, "d1", "df")
        self._breaking_labels = self._distance_labels_for(self.breaking_indices, "d0", "db")
        self._mode = self._infer_mode()
        if self.default_print_args is None:
            self.default_print_args = self._default_print_args_for_mode()

    @staticmethod
    def _distance_labels_for(
        pairs: Sequence[Tuple[int, int]],
        first_label: str,
        extra_prefix: str,
    ) -> List[str]:
        labels = []
        for i in range(len(pairs)):
            labels.append(first_label if i == 0 else f"{extra_prefix}{i}")
        return labels

    def _all_distance_labels(self) -> List[str]:
        return list(self._breaking_labels) + list(self._forming_labels)

    def _needs_sort(self) -> bool:
        n_distances = len(self.forming_indices) + len(self.breaking_indices)
        if self._mode == "difference":
            return True
        return n_distances >= 2

    def _wall_arg(self) -> str:
        n_distances = len(self.forming_indices) + len(self.breaking_indices)
        if self._mode == "difference":
            return "dsort.1"
        if n_distances >= 2:
            return f"dsort.{n_distances}"
        return self._all_distance_labels()[0]

    def _infer_mode(self) -> str:
        if self.forming_bonds and self.breaking_bonds:
            return "difference"
        if self.forming_bonds:
            return "forming"
        return "breaking"

    def _default_print_args_for_mode(self) -> str:
        """Build PRINT ARG list using the actual main CV label."""
        parts = self._all_distance_labels()
        if self._needs_sort():
            parts.append("dsort.1")
            wall_arg = self._wall_arg()
            if wall_arg not in parts:
                parts.append(wall_arg)
        cv_name = self.default_cv_name
        if cv_name not in parts:
            parts.append(cv_name)
        parts.append("mr.*")
        return ",".join(parts)

    def get_indices(self) -> Dict[str, int]:
        indices: Dict[str, int] = {}
        if self.index_forming_a is not None and self.index_forming_b is not None:
            indices["forming_a"] = self.index_forming_a
            indices["forming_b"] = self.index_forming_b
        if self.index_breaking_a is not None and self.index_breaking_b is not None:
            indices["breaking_a"] = self.index_breaking_a
            indices["breaking_b"] = self.index_breaking_b
        for i, (index_a, index_b) in enumerate(self.forming_indices):
            indices[f"forming_{i}_a"] = index_a
            indices[f"forming_{i}_b"] = index_b
        for i, (index_a, index_b) in enumerate(self.breaking_indices):
            indices[f"breaking_{i}_a"] = index_a
            indices[f"breaking_{i}_b"] = index_b
        return indices

    def _plumed_atoms(self, index_a: int, index_b: int) -> str:
        return (
            f"{to_plumed_index(index_a, self.idx_start_from)},"
            f"{to_plumed_index(index_b, self.idx_start_from)}"
        )

    def _ase_distance(self, index_a: int, index_b: int) -> float:
        return self.system.get_distance(
            to_ase_index(index_a, self.idx_start_from),
            to_ase_index(index_b, self.idx_start_from),
        )

    def _distance_lines(self) -> List[str]:
        lines = []
        for label, (index_a, index_b) in zip(self._breaking_labels, self.breaking_indices):
            lines.append(
                f"{label}: DISTANCE ATOMS={self._plumed_atoms(index_a, index_b)} NOPBC"
            )
        for label, (index_a, index_b) in zip(self._forming_labels, self.forming_indices):
            lines.append(
                f"{label}: DISTANCE ATOMS={self._plumed_atoms(index_a, index_b)} NOPBC"
            )
        return lines

    def _combine_line(self, cv_name: str) -> Optional[str]:
        args = list(self._forming_labels) + list(self._breaking_labels)
        if not args:
            return None
        if (
            len(args) == 1
            and self._mode != "difference"
            and cv_name == args[0]
        ):
            return None
        coefficients = ",".join(
            ["1"] * len(self._forming_labels) + ["-1"] * len(self._breaking_labels)
        )
        return f"{cv_name}: COMBINE ARG={','.join(args)} COEFFICIENTS={coefficients} PERIODIC=NO"

    def define_main_rc(self) -> Tuple[str, str]:
        cv_name = self.default_cv_name
        lines = self._distance_lines()
        if self._needs_sort():
            lines.append(f"dsort: SORT ARG={','.join(self._all_distance_labels())}")
        combine = self._combine_line(cv_name)
        if combine is not None:
            lines.append(combine)

        if combine is None:
            return self._all_distance_labels()[0], "\n".join(lines)
        return cv_name, "\n".join(lines)

    def define_additional_restraints(self) -> List[str]:
        if self.max_bond_length is None:
            return []
        return [
            f"uwall: UPPER_WALLS ARG={self._wall_arg()} AT={self.max_bond_length} "
            f"KAPPA={1000 * kcal / mol}"
        ]

    def calc_main_rc(self) -> float:
        forming_sum = sum(
            self._ase_distance(index_a, index_b)
            for index_a, index_b in self.forming_indices
        )
        breaking_sum = sum(
            self._ase_distance(index_a, index_b)
            for index_a, index_b in self.breaking_indices
        )
        if self._mode == "difference":
            return forming_sum - breaking_sum
        if self._mode == "forming":
            return forming_sum
        return breaking_sum
