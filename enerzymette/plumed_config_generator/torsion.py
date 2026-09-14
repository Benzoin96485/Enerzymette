"""Four-atom torsion (dihedral) PLUMED generator.

The main CV is a PLUMED ``TORSION`` over four atoms resolved in Python
(explicit indices or PDB selectors).  The emitted PLUMED input uses only
1-based atom numbers — never ``MOLINFO`` or named protein torsions.

Scan intervals are literal: ``lower_bound=-π`` and ``upper_bound=π`` is a
full turn of width ``2π``, not a zero-width periodic wrap.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

from ase import Atoms

from enerzymette.plumed_config_generator._engine import (
    PlumedConfigGenerator,
    ReactionCoordinate,
)
from enerzymette.plumed_config_generator.atom_selection import (
    TorsionSpec,
    coerce_atom_spec,
    coerce_torsion,
    resolve_torsion_indices,
    to_ase_index,
    to_plumed_index,
)

TorsionLike = Union[TorsionSpec, Mapping[str, Any], Sequence[Any]]

_ALLOWED_ANGLE_UNITS = {"rad", "deg"}


class TorsionConfigGenerator(PlumedConfigGenerator):
    """PLUMED generator for a four-atom dihedral scan / steered MD."""

    default_cv_name: str = "phi"
    default_print_args: Optional[str] = "phi,mr.*"
    default_dump_interval: int = 20

    def __init__(
        self,
        system: Atoms,
        *,
        atoms: Optional[TorsionLike] = None,
        atom1: Any = None,
        atom2: Any = None,
        atom3: Any = None,
        atom4: Any = None,
        angle_unit: str = "rad",
        cv_name: Optional[str] = None,
        **kwargs,
    ) -> None:
        super().__init__(system, **kwargs)
        unit = str(angle_unit).strip().lower()
        if unit in {"radian", "radians"}:
            unit = "rad"
        elif unit in {"degree", "degrees"}:
            unit = "deg"
        if unit not in _ALLOWED_ANGLE_UNITS:
            raise ValueError(
                f"angle_unit must be 'rad' or 'deg'; got {angle_unit!r}"
            )
        self.angle_unit = unit
        if cv_name is not None:
            self.default_cv_name = cv_name
        self.default_print_args = f"{self.default_cv_name},mr.*"

        named = (atom1, atom2, atom3, atom4)
        has_named = any(value is not None for value in named)
        if atoms is not None and has_named:
            raise ValueError("Provide atoms or atom1..atom4, not both")
        if atoms is not None:
            spec = coerce_torsion(atoms, label="atoms")
        elif has_named:
            if any(value is None for value in named):
                raise ValueError(
                    "atom1, atom2, atom3, and atom4 are all required when "
                    "atoms is omitted"
                )
            spec = TorsionSpec(
                atom1=coerce_atom_spec(atom1, label="atom1"),
                atom2=coerce_atom_spec(atom2, label="atom2"),
                atom3=coerce_atom_spec(atom3, label="atom3"),
                atom4=coerce_atom_spec(atom4, label="atom4"),
            )
            spec.validate(label="torsion")
        else:
            raise ValueError(
                "Provide atoms (list/mapping) or all of atom1, atom2, atom3, atom4"
            )
        assert spec is not None
        self.torsion = spec
        self.index_atom1, self.index_atom2, self.index_atom3, self.index_atom4 = (
            resolve_torsion_indices(
                spec,
                idx_start_from=self.idx_start_from,
                reference_pdb_file=self.reference_pdb,
                n_atoms=len(system),
                label="torsion",
            )
        )

    def _to_radians(self, value: float) -> float:
        if self.angle_unit == "deg":
            return math.radians(float(value))
        return float(value)

    def _resolved_user_bounds(
        self,
        lower_bound: Optional[float],
        upper_bound: Optional[float],
    ) -> Tuple[float, float]:
        if lower_bound is None:
            lower_bound = -180.0 if self.angle_unit == "deg" else -math.pi
        if upper_bound is None:
            upper_bound = 180.0 if self.angle_unit == "deg" else math.pi
        return lower_bound, upper_bound

    def build_reaction_coordinate(
        self,
        *,
        lower_bound: Optional[float] = None,
        upper_bound: Optional[float] = None,
        dump_interval: Optional[int] = None,
        kappa: Optional[float] = None,
        print_args: Optional[str] = None,
        **kwargs,
    ) -> ReactionCoordinate:
        user_lower, user_upper = self._resolved_user_bounds(lower_bound, upper_bound)
        return super().build_reaction_coordinate(
            lower_bound=self._to_radians(user_lower),
            upper_bound=self._to_radians(user_upper),
            dump_interval=(
                self.default_dump_interval if dump_interval is None else dump_interval
            ),
            kappa=kappa,
            print_args=print_args,
            **kwargs,
        )

    def choose_scan_endpoints(
        self,
        rc: ReactionCoordinate,
        *,
        target_value: Optional[float] = None,
        target_initial_value: Optional[float] = None,
    ) -> Tuple[float, float]:
        """Scan the configured interval literally; never wrap ±π to zero width.

        Without an explicit target, ``x0``/``x1`` are the configured bounds
        (radians).  With a target, start at the current principal-value torsion
        and go to the target value without taking the shortest periodic arc.
        """
        if target_value is not None:
            return rc.initial_value, float(target_value)
        if target_initial_value is not None:
            return rc.initial_value, float(target_initial_value)
        return rc.lower_bound, rc.upper_bound

    def get_indices(self) -> Dict[str, int]:
        return {
            "atom1": self.index_atom1,
            "atom2": self.index_atom2,
            "atom3": self.index_atom3,
            "atom4": self.index_atom4,
        }

    def define_main_rc(self) -> Tuple[str, str]:
        cv_name = self.default_cv_name
        plumed_atoms = ",".join(
            str(to_plumed_index(index, self.idx_start_from))
            for index in (
                self.index_atom1,
                self.index_atom2,
                self.index_atom3,
                self.index_atom4,
            )
        )
        return cv_name, f"{cv_name}: TORSION ATOMS={plumed_atoms} NOPBC"

    def calc_main_rc(self) -> float:
        return math.radians(
            float(
                self.system.get_dihedral(
                    to_ase_index(self.index_atom1, self.idx_start_from),
                    to_ase_index(self.index_atom2, self.idx_start_from),
                    to_ase_index(self.index_atom3, self.idx_start_from),
                    to_ase_index(self.index_atom4, self.idx_start_from),
                    mic=False,
                )
            )
        )

    def scan(
        self,
        *,
        target_value: float,
        lower_bound: Optional[float] = None,
        upper_bound: Optional[float] = None,
        dump_interval: Optional[int] = None,
        **kwargs,
    ) -> List[str]:
        user_lower, user_upper = self._resolved_user_bounds(lower_bound, upper_bound)
        return super().scan(
            target_value=target_value,
            lower_bound=user_lower,
            upper_bound=user_upper,
            dump_interval=(
                self.default_dump_interval if dump_interval is None else dump_interval
            ),
            **kwargs,
        )

    def standard_steered_md(
        self,
        *,
        integrate_config: Optional[dict] = None,
        lower_bound: Optional[float] = None,
        upper_bound: Optional[float] = None,
        dump_interval: Optional[int] = None,
        **kwargs,
    ) -> List[str]:
        user_lower, user_upper = self._resolved_user_bounds(lower_bound, upper_bound)
        return super().standard_steered_md(
            integrate_config=integrate_config,
            lower_bound=user_lower,
            upper_bound=user_upper,
            dump_interval=(
                self.default_dump_interval if dump_interval is None else dump_interval
            ),
            **kwargs,
        )

    def standard_restrained_md(
        self,
        *,
        integrate_config: Optional[dict] = None,
        lower_bound: Optional[float] = None,
        upper_bound: Optional[float] = None,
        dump_interval: Optional[int] = None,
        **kwargs,
    ) -> List[str]:
        user_lower, user_upper = self._resolved_user_bounds(lower_bound, upper_bound)
        return super().standard_restrained_md(
            integrate_config=integrate_config,
            lower_bound=user_lower,
            upper_bound=user_upper,
            dump_interval=(
                self.default_dump_interval if dump_interval is None else dump_interval
            ),
            **kwargs,
        )

    def naive_steered_md(
        self,
        *,
        integrate_config: Optional[dict] = None,
        lower_bound: Optional[float] = None,
        upper_bound: Optional[float] = None,
        dump_interval: Optional[int] = None,
        warmup_steps: int,
        **kwargs,
    ) -> List[str]:
        user_lower, user_upper = self._resolved_user_bounds(lower_bound, upper_bound)
        return super().naive_steered_md(
            integrate_config=integrate_config,
            lower_bound=user_lower,
            upper_bound=user_upper,
            dump_interval=(
                self.default_dump_interval if dump_interval is None else dump_interval
            ),
            warmup_steps=warmup_steps,
            **kwargs,
        )
