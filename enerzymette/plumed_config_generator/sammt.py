from typing import List, Tuple, Optional
from ase import Atoms
from ase.units import kJ, mol, fs, kcal

from enerzymette.plumed_config_generator._engine import (
    ReactionCoordinate,
    generate_steered_md,
    generate_scan_restraint,
)

_SAMMT_PRINT_ARGS = "d0,d1,dsort.1,dd,mr.*"


def get_sammt_index(
    idx_start_from: int,
    reference_pdb_file: str,
    substrate: str,
    nucleophile: str,
) -> Tuple[int, int, int]:
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
                elif resname == substrate:
                    if atomname == nucleophile:
                        index_nucleophile = atom_count
                atom_count += 1
    return index_sulphur + idx_start_from, index_methyl_carbon + idx_start_from, index_nucleophile + idx_start_from


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
        index_sulphur, index_methyl_carbon, index_nucleophile = get_sammt_index(
            idx_start_from, reference_pdb_file, substrate, nucleophile
        )
    if index_sulphur is None or index_methyl_carbon is None or index_nucleophile is None:
        raise ValueError("Index of sulphur, methyl carbon, and nucleophile must be provided")
    return index_sulphur, index_methyl_carbon, index_nucleophile


def _build_sammt_preamble(
    system: Atoms,
    idx_start_from: int,
    index_sulphur: int,
    index_methyl_carbon: int,
    index_nucleophile: int,
    max_bond_length: Optional[float],
) -> Tuple[List[str], float]:
    plumed_config = []
    plumed_config.append(f"UNITS LENGTH=A TIME={0.001 / fs} ENERGY={1 / kJ * mol}")
    plumed_config.append(
        f"d0: DISTANCE ATOMS={index_sulphur + 1 - idx_start_from},{index_methyl_carbon + 1 - idx_start_from} NOPBC"
    )
    plumed_config.append(
        f"d1: DISTANCE ATOMS={index_methyl_carbon + 1 - idx_start_from},{index_nucleophile + 1 - idx_start_from} NOPBC"
    )
    plumed_config.append("dsort: SORT ARG=d0,d1")
    if max_bond_length is not None:
        plumed_config.append(f"uwall: UPPER_WALLS ARG=dsort.1 AT={max_bond_length} KAPPA={1000 * kcal / mol}")
    plumed_config.append("dd: COMBINE ARG=d1,d0 COEFFICIENTS=1,-1 PERIODIC=NO")
    current_d0 = system.get_distance(
        index_sulphur - idx_start_from, index_methyl_carbon - idx_start_from
    )
    current_d1 = system.get_distance(
        index_methyl_carbon - idx_start_from, index_nucleophile - idx_start_from
    )
    current_dd = current_d1 - current_d0
    return plumed_config, current_dd


def get_sammt_reaction_coordinate(
    system: Atoms,
    idx_start_from: int,
    dump_interval: int,
    upper_bound: float,
    lower_bound: float,
    reference_pdb_file: Optional[str] = None,
    substrate: Optional[str] = None,
    nucleophile: Optional[str] = None,
    index_sulphur: Optional[int] = None,
    index_methyl_carbon: Optional[int] = None,
    index_nucleophile: Optional[int] = None,
    max_bond_length: float = 3.0,
    **kwargs,
) -> ReactionCoordinate:
    index_sulphur, index_methyl_carbon, index_nucleophile = _resolve_sammt_indices(
        idx_start_from,
        reference_pdb_file,
        substrate,
        nucleophile,
        index_sulphur,
        index_methyl_carbon,
        index_nucleophile,
    )
    preamble, current_dd = _build_sammt_preamble(
        system,
        idx_start_from,
        index_sulphur,
        index_methyl_carbon,
        index_nucleophile,
        max_bond_length,
    )
    return ReactionCoordinate(
        preamble=preamble,
        cv_name="dd",
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        initial_value=current_dd,
        dump_interval=dump_interval,
        print_args=_SAMMT_PRINT_ARGS,
    )


def get_sammt_basics(
    system: Atoms,
    idx_start_from: int,
    dump_interval: int,
    reference_pdb_file: Optional[str] = None,
    substrate: Optional[str] = None,
    nucleophile: Optional[str] = None,
    index_sulphur: Optional[int] = None,
    index_methyl_carbon: Optional[int] = None,
    index_nucleophile: Optional[int] = None,
    max_bond_length: float = 3.0,
    **kwargs,
) -> Tuple[List[str], float]:
    index_sulphur, index_methyl_carbon, index_nucleophile = _resolve_sammt_indices(
        idx_start_from,
        reference_pdb_file,
        substrate,
        nucleophile,
        index_sulphur,
        index_methyl_carbon,
        index_nucleophile,
    )
    return _build_sammt_preamble(
        system,
        idx_start_from,
        index_sulphur,
        index_methyl_carbon,
        index_nucleophile,
        max_bond_length,
    )


def get_sammt_config(
    system: Atoms,
    integrate_config: dict,
    idx_start_from: int,
    dump_interval: int,
    upper_bound: float,
    lower_bound: float,
    reference_pdb_file: Optional[str] = None,
    substrate: Optional[str] = None,
    nucleophile: Optional[str] = None,
    index_sulphur: Optional[int] = None,
    index_methyl_carbon: Optional[int] = None,
    index_nucleophile: Optional[int] = None,
    max_bond_length: float = 3.0,
    **kwargs,
) -> List[str]:
    rc = get_sammt_reaction_coordinate(
        system,
        idx_start_from,
        dump_interval,
        upper_bound,
        lower_bound,
        reference_pdb_file=reference_pdb_file,
        substrate=substrate,
        nucleophile=nucleophile,
        index_sulphur=index_sulphur,
        index_methyl_carbon=index_methyl_carbon,
        index_nucleophile=index_nucleophile,
        max_bond_length=max_bond_length,
        **kwargs,
    )
    return generate_steered_md(rc, integrate_config)


def get_sammt_scan_config(
    system: Atoms,
    integrate_config: dict,
    idx_start_from: int,
    target_value: float,
    dump_interval: int,
    upper_bound: float,
    lower_bound: float,
    reference_pdb_file: Optional[str] = None,
    substrate: Optional[str] = None,
    nucleophile: Optional[str] = None,
    index_sulphur: Optional[int] = None,
    index_methyl_carbon: Optional[int] = None,
    index_nucleophile: Optional[int] = None,
    max_bond_length: float = 3.0,
    **kwargs,
) -> List[str]:
    rc = get_sammt_reaction_coordinate(
        system,
        idx_start_from,
        dump_interval,
        upper_bound,
        lower_bound,
        reference_pdb_file=reference_pdb_file,
        substrate=substrate,
        nucleophile=nucleophile,
        index_sulphur=index_sulphur,
        index_methyl_carbon=index_methyl_carbon,
        index_nucleophile=index_nucleophile,
        max_bond_length=None,
        **kwargs,
    )
    return generate_scan_restraint(rc, target_value)


def get_naive_sammt_config(
    system: Atoms,
    integrate_config: dict,
    idx_start_from: int,
    dump_interval: int,
    upper_bound: float,
    lower_bound: float,
    warmup_steps: int,
    reference_pdb_file: Optional[str] = None,
    substrate: Optional[str] = None,
    nucleophile: Optional[str] = None,
    index_sulphur: Optional[int] = None,
    index_methyl_carbon: Optional[int] = None,
    index_nucleophile: Optional[int] = None,
    **kwargs,
) -> List[str]:
    plumed_config, current_dd = get_sammt_basics(
        system,
        idx_start_from,
        dump_interval,
        reference_pdb_file,
        substrate,
        nucleophile,
        index_sulphur,
        index_methyl_carbon,
        index_nucleophile,
    )
    n_step = integrate_config.get("n_step")
    moving_restraint_kappa = 1000 * kcal / mol
    moving_restraint_component = (
        f"mr: MOVINGRESTRAINT ARG=dd STEP0=0 AT0={current_dd} KAPPA0={moving_restraint_kappa} "
        f"STEP1={warmup_steps} AT1={upper_bound} STEP2={n_step} AT2={lower_bound}"
    )
    plumed_config.append(moving_restraint_component)
    plumed_config.append(f"PRINT ARG=d0,d1,dsort.1,dd,mr.* STRIDE={dump_interval}")
    plumed_config.append(f"FLUSH STRIDE={dump_interval}")
    return plumed_config
