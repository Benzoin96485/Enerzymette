import json
import os
from typing import List, Tuple, Optional

from ase import Atoms
from ase.units import kJ, mol, fs, kcal
from mendeleev import element

from enerzymette.plumed_config_generator._engine import (
    ProtonTransferBias,
    ProtonTransferScope,
    ReactionCoordinate,
    append_proton_transfer_to_plumed,
    default_opes_state_wstride,
    generate_scan_restraint,
    generate_steered_md,
)

_SAMMT_PRINT_ARGS = "d0,d1,dsort.1,dd,mr.*"

_DEFAULT_PROTON_HBOND_COEFF = 1.1
_DEFAULT_ACCEPTOR_RADIUS = 4.0
_DEFAULT_OPES_PACE = 20
_DEFAULT_OPES_BARRIER_KJ_MOL = 20.0
_DEFAULT_OPES_BIASFACTOR = 100.0


def _covalent_radius_angstrom(symbol: str) -> float:
    """Return mendeleev covalent radius in Angstrom (mendeleev uses pm)."""
    radius_pm = element(symbol).covalent_radius
    if radius_pm is None:
        raise ValueError(f"No covalent radius available for element {symbol!r}")
    return radius_pm / 100.0


def _nearest_heavy_atom_index(system: Atoms, atom_idx: int) -> int:
    heavy_atoms = [i for i, atom in enumerate(system) if atom.number > 1]
    if not heavy_atoms:
        raise ValueError("No heavy atoms found in system")
    return min(
        heavy_atoms,
        key=lambda i: system.get_distance(atom_idx, i, mic=False),
    )


def _heavy_bond_counts_from_topology(
    topology_mol_file: str,
    expected_n_atoms: int,
) -> List[int]:
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ImportError(
            "RDKit is required to filter proton acceptors by topology_mol_file"
        ) from exc

    mol = Chem.MolFromMolFile(topology_mol_file, removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError(f"Could not read topology mol file: {topology_mol_file}")
    if mol.GetNumAtoms() != expected_n_atoms:
        raise ValueError(
            f"Topology atom count ({mol.GetNumAtoms()}) does not match "
            f"XYZ atom count ({expected_n_atoms}) for {topology_mol_file}"
        )

    heavy_bond_counts = [0 for _ in range(expected_n_atoms)]
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum() <= 1:
            continue
        heavy_bond_counts[idx] = sum(
            1
            for bond in atom.GetBonds()
            if bond.GetOtherAtom(atom).GetAtomicNum() > 1
        )
    return heavy_bond_counts


def _is_topologically_allowed_acceptor(
    atomic_number: int,
    heavy_bond_count: int,
) -> bool:
    if atomic_number == 8:
        return heavy_bond_count < 2
    if atomic_number == 7:
        return heavy_bond_count < 3
    return False


def resolve_proton_transfer_scope(
    system: Atoms,
    index_nucleophile: int,
    idx_start_from: int,
    proton_hbond_coeff: float = _DEFAULT_PROTON_HBOND_COEFF,
    acceptor_radius: float = _DEFAULT_ACCEPTOR_RADIUS,
    topology_mol_file: Optional[str] = None,
) -> ProtonTransferScope:
    """Identify proton donor, transfer protons, and acceptors for SAM-MT."""
    donor_ase_idx = index_nucleophile - idx_start_from
    if donor_ase_idx < 0 or donor_ase_idx >= len(system):
        raise ValueError(
            f"index_nucleophile={index_nucleophile} is out of range for {len(system)} atoms"
        )

    donor_symbol = system[donor_ase_idx].symbol
    proton_threshold = (
        _covalent_radius_angstrom(donor_symbol) + _covalent_radius_angstrom("H")
    ) * proton_hbond_coeff

    transfer_protons = [
        i
        for i, atom in enumerate(system)
        if atom.number == 1
        and system.get_distance(donor_ase_idx, i, mic=False) < proton_threshold
        and _nearest_heavy_atom_index(system, i) == donor_ase_idx
    ]
    if not transfer_protons:
        raise ValueError(
            f"No transfer protons found within {proton_threshold:.3f} A of donor "
            f"atom {donor_ase_idx} ({donor_symbol}) whose nearest heavy atom is "
            "the donor"
        )

    heavy_bond_counts = (
        _heavy_bond_counts_from_topology(topology_mol_file, len(system))
        if topology_mol_file is not None
        else None
    )
    acceptor_candidates = [
        i
        for i, atom in enumerate(system)
        if atom.number in (7, 8) and i != donor_ase_idx
        and (
            heavy_bond_counts is None
            or _is_topologically_allowed_acceptor(
                atom.number,
                heavy_bond_counts[i],
            )
        )
    ]
    if not acceptor_candidates:
        raise ValueError("No N/O proton acceptor candidates found near the donor")

    acceptors = [
        i
        for i in acceptor_candidates
        if system.get_distance(donor_ase_idx, i, mic=False) < acceptor_radius
    ]
    if not acceptors:
        raise ValueError(
            f"No acceptors found within {acceptor_radius:.3f} A of donor "
            f"atom {donor_ase_idx}"
        )

    return ProtonTransferScope(
        index_donor=donor_ase_idx,
        transfer_protons=transfer_protons,
        acceptors=acceptors,
    )


def load_proton_transfer_scope(scope_file: str) -> ProtonTransferScope:
    with open(scope_file, "r") as handle:
        return ProtonTransferScope.from_dict(json.load(handle))


def save_proton_transfer_scope(scope_file: str, scope: ProtonTransferScope) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(scope_file)), exist_ok=True)
    with open(scope_file, "w") as handle:
        json.dump(scope.to_dict(), handle, indent=2)


def _load_or_resolve_proton_transfer_scope(
    system: Atoms,
    index_nucleophile: int,
    idx_start_from: int,
    pt_scope_file: Optional[str],
    proton_hbond_coeff: float,
    acceptor_radius: float,
    topology_mol_file: Optional[str],
) -> ProtonTransferScope:
    if pt_scope_file is not None and os.path.exists(pt_scope_file):
        return load_proton_transfer_scope(pt_scope_file)

    scope = resolve_proton_transfer_scope(
        system,
        index_nucleophile,
        idx_start_from,
        proton_hbond_coeff=proton_hbond_coeff,
        acceptor_radius=acceptor_radius,
        topology_mol_file=topology_mol_file,
    )
    if pt_scope_file is not None:
        save_proton_transfer_scope(pt_scope_file, scope)
    return scope


def _maybe_append_proton_transfer(
    plumed_config: List[str],
    system: Atoms,
    idx_start_from: int,
    index_nucleophile: int,
    dump_interval: int,
    integrate_config: dict,
    proton_transfer: bool,
    pt_scope_file: Optional[str],
    pt_restart: bool,
    pt_state_file: str,
    proton_hbond_coeff: float,
    acceptor_radius: float,
    topology_mol_file: Optional[str],
    opes_pace: int,
    opes_barrier: float,
    opes_biasfactor: float,
    opes_state_wstride: Optional[int] = None,
) -> List[str]:
    if not proton_transfer:
        return plumed_config

    scope = _load_or_resolve_proton_transfer_scope(
        system,
        index_nucleophile,
        idx_start_from,
        pt_scope_file,
        proton_hbond_coeff,
        acceptor_radius,
        topology_mol_file,
    )
    n_step = integrate_config.get("n_step")
    state_wstride = (
        opes_state_wstride
        if opes_state_wstride is not None
        else default_opes_state_wstride(n_step, dump_interval)
    )
    bias = ProtonTransferBias(
        pace=opes_pace,
        barrier_kj_mol=opes_barrier,
        biasfactor=opes_biasfactor,
        state_wfile=pt_state_file,
        state_wstride=state_wstride,
        state_rfile=pt_state_file if pt_restart else None,
        pt_restart=pt_restart,
    )
    return append_proton_transfer_to_plumed(
        plumed_config,
        scope,
        idx_start_from,
        dump_interval,
        bias=bias,
    )


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
    proton_transfer: bool = False,
    pt_scope_file: Optional[str] = None,
    pt_restart: bool = False,
    pt_state_file: str = "opes_state.data",
    proton_hbond_coeff: float = _DEFAULT_PROTON_HBOND_COEFF,
    acceptor_radius: float = _DEFAULT_ACCEPTOR_RADIUS,
    topology_mol_file: Optional[str] = None,
    opes_pace: int = _DEFAULT_OPES_PACE,
    opes_barrier: float = _DEFAULT_OPES_BARRIER_KJ_MOL,
    opes_biasfactor: float = _DEFAULT_OPES_BIASFACTOR,
    opes_state_wstride: Optional[int] = None,
    **kwargs,
) -> List[str]:
    index_sulphur, index_methyl_carbon, index_nucleophile = _resolve_sammt_indices(
        idx_start_from,
        reference_pdb_file,
        substrate,
        nucleophile,
        index_sulphur,
        index_methyl_carbon,
        index_nucleophile,
    )
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
    )
    plumed_config = generate_steered_md(rc, integrate_config)
    return _maybe_append_proton_transfer(
        plumed_config,
        system,
        idx_start_from,
        index_nucleophile,
        dump_interval,
        integrate_config,
        proton_transfer,
        pt_scope_file,
        pt_restart,
        pt_state_file,
        proton_hbond_coeff,
        acceptor_radius,
        topology_mol_file,
        opes_pace,
        opes_barrier,
        opes_biasfactor,
        opes_state_wstride,
    )


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
