"""Generic forming / breaking bond reaction-coordinate generator.

Users supply up to two complete atom pairs:

* ``forming_bond``: the bond being formed
* ``breaking_bond``: the bond being broken

The main collective variable is chosen automatically:

* forming only  → ``rc = d_forming``
* breaking only → ``rc = d_breaking``
* both         → ``rc = d_forming - d_breaking``

Atoms may be given as explicit indices or as PDB selectors; see
:mod:`enerzymette.plumed_config_generator.atom_selection`.
"""

from __future__ import annotations

from typing import Any, Dict, Mapping, Optional, Tuple, Union

from ase import Atoms
from ase.units import kcal, mol

from enerzymette.plumed_config_generator._engine import PlumedConfigGenerator
from enerzymette.plumed_config_generator.atom_selection import (
    BondPairSpec,
    coerce_bond_pair,
    resolve_bond_pair_indices,
    to_ase_index,
    to_plumed_index,
)

BondPairInput = Optional[Union[BondPairSpec, Mapping[str, Any]]]


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
        max_bond_length: Optional[float] = None,
        cv_name: Optional[str] = None,
        **kwargs,
    ) -> None:
        super().__init__(system, **kwargs)
        self.forming_bond = coerce_bond_pair(forming_bond, label="forming_bond")
        self.breaking_bond = coerce_bond_pair(breaking_bond, label="breaking_bond")
        if self.forming_bond is None and self.breaking_bond is None:
            raise ValueError(
                "At least one of forming_bond or breaking_bond must be provided"
            )
        self.max_bond_length = max_bond_length
        if cv_name is not None:
            self.default_cv_name = cv_name

        n_atoms = len(system)
        pdb_file = self.reference_pdb
        self.index_forming_a: Optional[int] = None
        self.index_forming_b: Optional[int] = None
        self.index_breaking_a: Optional[int] = None
        self.index_breaking_b: Optional[int] = None

        if self.forming_bond is not None:
            self.index_forming_a, self.index_forming_b = resolve_bond_pair_indices(
                self.forming_bond,
                idx_start_from=self.idx_start_from,
                reference_pdb_file=pdb_file,
                n_atoms=n_atoms,
                label="forming_bond",
            )
        if self.breaking_bond is not None:
            self.index_breaking_a, self.index_breaking_b = resolve_bond_pair_indices(
                self.breaking_bond,
                idx_start_from=self.idx_start_from,
                reference_pdb_file=pdb_file,
                n_atoms=n_atoms,
                label="breaking_bond",
            )

        self._mode = self._infer_mode()
        if self.default_print_args is None:
            self.default_print_args = self._default_print_args_for_mode(self._mode)

    def _infer_mode(self) -> str:
        if self.forming_bond is not None and self.breaking_bond is not None:
            return "difference"
        if self.forming_bond is not None:
            return "forming"
        return "breaking"

    @staticmethod
    def _default_print_args_for_mode(mode: str) -> str:
        if mode == "difference":
            return "d0,d1,dsort.1,rc,mr.*"
        if mode == "forming":
            return "d1,rc,mr.*"
        return "d0,rc,mr.*"

    def get_indices(self) -> Dict[str, int]:
        indices: Dict[str, int] = {}
        if self.index_forming_a is not None and self.index_forming_b is not None:
            indices["forming_a"] = self.index_forming_a
            indices["forming_b"] = self.index_forming_b
        if self.index_breaking_a is not None and self.index_breaking_b is not None:
            indices["breaking_a"] = self.index_breaking_a
            indices["breaking_b"] = self.index_breaking_b
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

    def define_main_rc(self) -> Tuple[str, str]:
        cv_name = self.default_cv_name
        lines = []

        if self._mode in {"breaking", "difference"}:
            assert self.index_breaking_a is not None and self.index_breaking_b is not None
            lines.append(
                f"d0: DISTANCE ATOMS={self._plumed_atoms(self.index_breaking_a, self.index_breaking_b)} NOPBC"
            )
        if self._mode in {"forming", "difference"}:
            assert self.index_forming_a is not None and self.index_forming_b is not None
            lines.append(
                f"d1: DISTANCE ATOMS={self._plumed_atoms(self.index_forming_a, self.index_forming_b)} NOPBC"
            )

        if self._mode == "difference":
            lines.append("dsort: SORT ARG=d0,d1")
            if self.max_bond_length is not None:
                lines.append(
                    f"uwall: UPPER_WALLS ARG=dsort.1 AT={self.max_bond_length} "
                    f"KAPPA={1000 * kcal / mol}"
                )
            lines.append(f"{cv_name}: COMBINE ARG=d1,d0 COEFFICIENTS=1,-1 PERIODIC=NO")
        elif self._mode == "forming":
            if self.max_bond_length is not None:
                lines.append(
                    f"uwall: UPPER_WALLS ARG=d1 AT={self.max_bond_length} "
                    f"KAPPA={1000 * kcal / mol}"
                )
            if cv_name != "d1":
                lines.append(f"{cv_name}: COMBINE ARG=d1 COEFFICIENTS=1 PERIODIC=NO")
            else:
                # CV is already named d1.
                pass
        else:  # breaking
            if self.max_bond_length is not None:
                lines.append(
                    f"uwall: UPPER_WALLS ARG=d0 AT={self.max_bond_length} "
                    f"KAPPA={1000 * kcal / mol}"
                )
            if cv_name != "d0":
                lines.append(f"{cv_name}: COMBINE ARG=d0 COEFFICIENTS=1 PERIODIC=NO")

        # For single-bond modes when cv_name equals the distance label, ensure
        # define_main_rc still returns a usable name without a redundant COMBINE.
        if self._mode == "forming" and cv_name == "d1":
            return "d1", "\n".join(lines)
        if self._mode == "breaking" and cv_name == "d0":
            return "d0", "\n".join(lines)
        return cv_name, "\n".join(lines)

    def calc_main_rc(self) -> float:
        if self._mode == "difference":
            assert self.index_breaking_a is not None and self.index_breaking_b is not None
            assert self.index_forming_a is not None and self.index_forming_b is not None
            d0 = self._ase_distance(self.index_breaking_a, self.index_breaking_b)
            d1 = self._ase_distance(self.index_forming_a, self.index_forming_b)
            return d1 - d0
        if self._mode == "forming":
            assert self.index_forming_a is not None and self.index_forming_b is not None
            return self._ase_distance(self.index_forming_a, self.index_forming_b)
        assert self.index_breaking_a is not None and self.index_breaking_b is not None
        return self._ase_distance(self.index_breaking_a, self.index_breaking_b)
