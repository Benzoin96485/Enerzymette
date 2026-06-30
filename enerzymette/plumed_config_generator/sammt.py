from typing import Dict, Optional, Tuple

from ase import Atoms
from ase.units import kcal, mol

from enerzymette.plumed_config_generator._engine import PlumedConfigGenerator

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
    """Resolve SAM SD, SAM CE, and substrate nucleophile indices from a PDB file."""
    index_sulphur = None
    index_methyl_carbon = None
    index_nucleophile = None

    with open(reference_pdb_file, "r") as f:
        atom_count = 0
        for line in f.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resname = line[17:20].strip()
                atomname = line[11:16].strip()
                if resname == "SAM":
                    if atomname == "SD":
                        index_sulphur = atom_count
                    elif atomname == "CE":
                        index_methyl_carbon = atom_count
                elif resname == substrate and atomname == nucleophile:
                    index_nucleophile = atom_count
                atom_count += 1

    if index_sulphur is None:
        raise ValueError(f"Could not find SAM SD in {reference_pdb_file}")
    if index_methyl_carbon is None:
        raise ValueError(f"Could not find SAM CE in {reference_pdb_file}")
    if index_nucleophile is None:
        raise ValueError(
            f"Could not find nucleophile atom {nucleophile!r} in residue "
            f"{substrate!r} from {reference_pdb_file}"
        )

    return (
        index_sulphur + idx_start_from,
        index_methyl_carbon + idx_start_from,
        index_nucleophile + idx_start_from,
    )


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


class SAMMTConfigGenerator(PlumedConfigGenerator):
    """PLUMED generator for SAM-dependent methyltransferase reactions."""

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
        **kwargs,
    ) -> None:
        super().__init__(system, **kwargs)
        self.substrate = substrate
        self.nucleophile = nucleophile
        self.max_bond_length = max_bond_length
        (
            self.index_sulphur,
            self.index_methyl_carbon,
            self.index_nucleophile,
        ) = _resolve_sammt_indices(
            self.idx_start_from,
            self.reference_pdb,
            substrate,
            nucleophile,
            index_sulphur,
            index_methyl_carbon,
            index_nucleophile,
        )
        if self.proton_transfer_config.enabled and self.proton_transfer_config.donor is None:
            self.proton_transfer_config.donor = "nucleophile"

    def get_indices(self) -> Dict[str, int]:
        return {
            "sulphur": self.index_sulphur,
            "sulfur": self.index_sulphur,
            "methyl_carbon": self.index_methyl_carbon,
            "nucleophile": self.index_nucleophile,
        }

    def define_main_rc(self) -> Tuple[str, str]:
        lines = [
            (
                f"d0: DISTANCE ATOMS={self.index_sulphur + 1 - self.idx_start_from},"
                f"{self.index_methyl_carbon + 1 - self.idx_start_from} NOPBC"
            ),
            (
                f"d1: DISTANCE ATOMS={self.index_methyl_carbon + 1 - self.idx_start_from},"
                f"{self.index_nucleophile + 1 - self.idx_start_from} NOPBC"
            ),
            "dsort: SORT ARG=d0,d1",
        ]
        if self.max_bond_length is not None:
            lines.append(
                f"uwall: UPPER_WALLS ARG=dsort.1 AT={self.max_bond_length} "
                f"KAPPA={1000 * kcal / mol}"
            )
        lines.append("dd: COMBINE ARG=d1,d0 COEFFICIENTS=1,-1 PERIODIC=NO")
        return "dd", "\n".join(lines)

    def calc_main_rc(self) -> float:
        current_d0 = self.system.get_distance(
            self.index_sulphur - self.idx_start_from,
            self.index_methyl_carbon - self.idx_start_from,
        )
        current_d1 = self.system.get_distance(
            self.index_methyl_carbon - self.idx_start_from,
            self.index_nucleophile - self.idx_start_from,
        )
        return current_d1 - current_d0
