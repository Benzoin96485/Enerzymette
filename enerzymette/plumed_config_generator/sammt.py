"""SAM-dependent methyltransferase (SAMMT) PLUMED generator.

This is a specialized :class:`BondReactionConfigGenerator` with:

* breaking bond = SAM SD – SAM CE
* forming bond  = SAM CE – substrate nucleophile
* main CV ``dd = d(CE, Nu) - d(SD, CE)``

Users typically provide only ``substrate`` and ``nucleophile`` together with
``reference_pdb_file``; explicit ``index_sulphur`` / ``index_methyl_carbon`` /
``index_nucleophile`` remain supported.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

from ase import Atoms

from enerzymette.plumed_config_generator.atom_selection import (
    AtomSpec,
    BondPairSpec,
    parse_pdb_atoms,
    resolve_atom_index,
)
from enerzymette.plumed_config_generator.bond_reaction import BondReactionConfigGenerator

_SAMMT_PRINT_ARGS = "d0,d1,dsort.1,dd,mr.*"


def get_sammt_scan_bond_indices(
    reference_pdb_file: str,
    substrate: str,
    nucleophile: str,
) -> Tuple[int, int]:
    """Return 0-based methyl-carbon and nucleophile indices for ASE bond scans."""
    _, index_methyl_carbon, index_nucleophile = get_sammt_index(
        0,
        reference_pdb_file,
        substrate,
        nucleophile,
    )
    return index_methyl_carbon, index_nucleophile


def get_sammt_index(
    idx_start_from: int,
    reference_pdb_file: str,
    substrate: str,
    nucleophile: str,
) -> Tuple[int, int, int]:
    """Resolve SAM SD, SAM CE, and substrate nucleophile indices from a PDB file.

    Matches the historical SAMMT convention: locate the unique SAM residue
    atoms ``SD`` / ``CE`` and the unique ``(substrate, nucleophile)`` atom by
    residue name and atom name (chain ID and residue number are not required
    when the match is unique).
    """
    records = parse_pdb_atoms(reference_pdb_file)
    index_sulphur = resolve_atom_index(
        AtomSpec(residue_name="SAM", atom_name="SD"),
        idx_start_from=idx_start_from,
        reference_pdb_file=reference_pdb_file,
        pdb_records=records,
        label="SAM SD",
    )
    index_methyl_carbon = resolve_atom_index(
        AtomSpec(residue_name="SAM", atom_name="CE"),
        idx_start_from=idx_start_from,
        reference_pdb_file=reference_pdb_file,
        pdb_records=records,
        label="SAM CE",
    )
    index_nucleophile = resolve_atom_index(
        AtomSpec(residue_name=substrate, atom_name=nucleophile),
        idx_start_from=idx_start_from,
        reference_pdb_file=reference_pdb_file,
        pdb_records=records,
        label=f"nucleophile {nucleophile!r} in residue {substrate!r}",
    )
    return index_sulphur, index_methyl_carbon, index_nucleophile


def _resolve_sammt_indices(
    idx_start_from: int,
    reference_pdb_file: Optional[str],
    substrate: Optional[str],
    nucleophile: Optional[str],
    index_sulphur: Optional[int],
    index_methyl_carbon: Optional[int],
    index_nucleophile: Optional[int],
) -> Tuple[int, int, int]:
    if reference_pdb_file is not None:
        if substrate is None or nucleophile is None:
            raise ValueError(
                "substrate and nucleophile are required with reference_pdb_file"
            )
        index_sulphur, index_methyl_carbon, index_nucleophile = get_sammt_index(
            idx_start_from, reference_pdb_file, substrate, nucleophile
        )
    if index_sulphur is None or index_methyl_carbon is None or index_nucleophile is None:
        raise ValueError(
            "index_sulphur, index_methyl_carbon, and index_nucleophile must be "
            "provided when reference_pdb_file is omitted"
        )
    return index_sulphur, index_methyl_carbon, index_nucleophile


class SAMMTConfigGenerator(BondReactionConfigGenerator):
    """PLUMED generator for SAM-dependent methyltransferase reactions."""

    default_cv_name = "dd"
    default_print_args = _SAMMT_PRINT_ARGS

    def __init__(
        self,
        system: Atoms,
        *,
        substrate: Optional[str] = None,
        nucleophile: Optional[str] = None,
        index_sulphur: Optional[int] = None,
        index_methyl_carbon: Optional[int] = None,
        index_nucleophile: Optional[int] = None,
        max_bond_length: Optional[float] = 3.0,
        forming_bond=None,
        breaking_bond=None,
        **kwargs,
    ) -> None:
        if forming_bond is not None or breaking_bond is not None:
            raise ValueError(
                "SAMMTConfigGenerator does not accept forming_bond/breaking_bond; "
                "use substrate/nucleophile or explicit SAMMT indices instead"
            )

        idx_start_from = int(kwargs.get("idx_start_from", 1))
        reference_pdb = kwargs.get("reference_pdb_file") or kwargs.get("reference_pdb")
        (
            index_sulphur,
            index_methyl_carbon,
            index_nucleophile,
        ) = _resolve_sammt_indices(
            idx_start_from,
            reference_pdb,
            substrate,
            nucleophile,
            index_sulphur,
            index_methyl_carbon,
            index_nucleophile,
        )

        super().__init__(
            system,
            forming_bond=BondPairSpec(
                atom1=AtomSpec(index=index_methyl_carbon),
                atom2=AtomSpec(index=index_nucleophile),
            ),
            breaking_bond=BondPairSpec(
                atom1=AtomSpec(index=index_sulphur),
                atom2=AtomSpec(index=index_methyl_carbon),
            ),
            max_bond_length=max_bond_length,
            **kwargs,
        )
        self.substrate = substrate
        self.nucleophile = nucleophile
        self.index_sulphur = index_sulphur
        self.index_methyl_carbon = index_methyl_carbon
        self.index_nucleophile = index_nucleophile

        if self.proton_transfer_config.enabled and self.proton_transfer_config.donor is None:
            self.proton_transfer_config.donor = "nucleophile"

    def get_indices(self) -> Dict[str, int]:
        return {
            "sulphur": self.index_sulphur,
            "sulfur": self.index_sulphur,
            "methyl_carbon": self.index_methyl_carbon,
            "nucleophile": self.index_nucleophile,
            "forming_a": self.index_methyl_carbon,
            "forming_b": self.index_nucleophile,
            "breaking_a": self.index_sulphur,
            "breaking_b": self.index_methyl_carbon,
        }
