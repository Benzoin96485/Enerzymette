from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Optional, Tuple, Union
import json
import warnings

from ase import Atoms
from ase.units import kJ, kcal, mol

if TYPE_CHECKING:
    from ._engine import PlumedConfigGenerator


AtomSelector = Union[str, int]


@dataclass
class ProtonTransferConfig:
    """Configuration for a proton-transfer PLUMED plugin.

    Atom selectors may be descriptive names returned by
    :meth:`PlumedConfigGenerator.get_indices` or explicit atom indices in the
    generator's ``idx_start_from`` convention.
    """

    enabled: bool = False
    plugin: str = "local_opes"
    donor: Optional[AtomSelector] = None
    proton: Optional[Union[AtomSelector, List[AtomSelector]]] = None
    acceptor: Optional[Union[AtomSelector, List[AtomSelector]]] = None
    flavor: Optional[str] = None
    scope_file: Optional[str] = None
    restart: bool = False
    state_file: str = "opes_state.data"
    proton_hbond_coeff: float = 1.1
    acceptor_radius: float = 4.0
    topology_mol_file: Optional[str] = None
    opes_pace: int = 20
    opes_barrier: float = 20.0
    opes_biasfactor: float = 100.0
    opes_state_wstride: Optional[int] = None
    acceptor_grouping: str = "element"
    coordination_r0: float = 2.0
    coordination_d0: float = 0.1
    coupled_pt_weight: float = 1.0
    coupled_env_weight: float = 0.25
    coupled_dd_weight: float = 0.1
    geometry_walls: bool = True
    max_donor_h_distance: float = 2.5
    max_donor_acceptor_distance: float = 6.0
    single_proton: bool = True


class ProtonTransferPlugin(ABC):
    """Base class for pluggable proton-transfer CV/bias builders."""

    name: str

    @abstractmethod
    def append_to_plumed(
        self,
        generator: "PlumedConfigGenerator",
        plumed_config: List[str],
        config: ProtonTransferConfig,
        dump_interval: int,
        integrate_config: dict,
    ) -> List[str]:
        """Return a PLUMED config with this plugin inserted."""


PROTON_TRANSFER_PLUGINS: Dict[str, ProtonTransferPlugin] = {}


def register_proton_transfer_plugin(plugin: ProtonTransferPlugin) -> None:
    PROTON_TRANSFER_PLUGINS[plugin.name] = plugin


def build_proton_transfer_config(
    proton_transfer: Optional[Union[bool, dict, ProtonTransferConfig]],
    *,
    proton_transfer_plugin: Optional[str] = None,
    topology_mol_file: Optional[str] = None,
    **kwargs,
) -> ProtonTransferConfig:
    if isinstance(proton_transfer, ProtonTransferConfig):
        config = proton_transfer
        if proton_transfer_plugin is not None:
            config.plugin = proton_transfer_plugin
        return config

    if isinstance(proton_transfer, dict):
        data = dict(proton_transfer)
        data.setdefault("enabled", True)
    else:
        data = {"enabled": bool(proton_transfer)}

    if topology_mol_file is not None and "topology_mol_file" not in data:
        data["topology_mol_file"] = topology_mol_file
    if proton_transfer_plugin is not None:
        data["plugin"] = proton_transfer_plugin
    return ProtonTransferConfig(**data)


def append_optional_proton_transfer(
    generator: "PlumedConfigGenerator",
    plumed_config: List[str],
    *,
    dump_interval: int,
    integrate_config: Optional[dict] = None,
    proton_transfer: Optional[Union[bool, dict, ProtonTransferConfig]] = None,
    proton_transfer_plugin: Optional[str] = None,
) -> List[str]:
    config = generator.proton_transfer_config
    if proton_transfer is not None or proton_transfer_plugin is not None:
        config = build_proton_transfer_config(
            proton_transfer if proton_transfer is not None else config,
            proton_transfer_plugin=proton_transfer_plugin,
        )
    if not config.enabled:
        return plumed_config

    plugin = PROTON_TRANSFER_PLUGINS.get(config.plugin)
    if plugin is None:
        warnings.warn(
            f"Unknown proton-transfer plugin {config.plugin!r}; skipping.",
            RuntimeWarning,
        )
        return plumed_config
    try:
        return plugin.append_to_plumed(
            generator,
            plumed_config,
            config,
            dump_interval,
            integrate_config or generator.integrate_config,
        )
    except Exception as exc:
        warnings.warn(
            f"Could not insert proton-transfer plugin {config.plugin!r}: {exc}",
            RuntimeWarning,
        )
        return plumed_config


def _as_list(value):
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        return list(value)
    return [value]


@dataclass
class ProtonTransferScope:
    """ASE 0-based atom indices for proton-transfer CV."""

    index_donor: int
    transfer_protons: List[int]
    acceptors: List[int]

    def to_dict(self) -> dict:
        return {
            "index_donor": self.index_donor,
            "transfer_protons": list(self.transfer_protons),
            "acceptors": list(self.acceptors),
        }

    @classmethod
    def from_dict(cls, data: dict) -> "ProtonTransferScope":
        return cls(
            index_donor=int(data["index_donor"]),
            transfer_protons=[int(i) for i in data["transfer_protons"]],
            acceptors=[int(i) for i in data["acceptors"]],
        )


def default_opes_state_wstride(n_step: Optional[int], dump_interval: int) -> int:
    """Pick a STATE_WSTRIDE that divides n_step so the final MD step is checkpointed.

    ASE/Enerzyme do not emit PLUMED CPT events; without STATE_WSTRIDE, STATE_WFILE
    is never written and OPES restart cannot work.
    """
    if n_step is None or n_step <= 0:
        return 500
    for candidate in (1000, 500, 200, 100, dump_interval, n_step):
        if 0 < candidate <= n_step and n_step % candidate == 0:
            return candidate
    return n_step


@dataclass
class ProtonTransferBias:
    pace: int = 20
    barrier_kj_mol: float = 1.0
    biasfactor: float = 15.0
    state_wfile: str = "opes_state.data"
    state_wstride: int = 500
    state_rfile: Optional[str] = None
    restart: bool = False
    min_beta: float = 50.0


@dataclass
class ProtonTransferCVConfig:
    """Generic local proton-network CV options (SAM-MT transferable)."""

    mode: str = "nearest_distance"
    acceptor_grouping: str = "element"
    coordination_r0: float = 2.0
    coordination_d0: float = 0.1
    coupled_pt_weight: float = 1.0
    coupled_env_weight: float = 0.25
    coupled_dd_weight: float = 0.1
    geometry_walls: bool = True
    max_donor_h_distance: float = 2.5
    max_donor_acceptor_distance: float = 6.0
    single_proton: bool = True


def _plumed_atom_list(indices: List[int], idx_start_from: int) -> str:
    return ",".join(str(i + 1 - idx_start_from) for i in indices)


def _resolve_scope_atoms(
    scope: ProtonTransferScope,
    cv_config: Optional[ProtonTransferCVConfig],
    system: Optional[Atoms],
) -> Tuple[List[int], List[int]]:
    protons = list(scope.transfer_protons)
    acceptors = list(scope.acceptors)
    if cv_config.single_proton and len(protons) > 1 and system is not None:
        donor = scope.index_donor
        protons = [
            min(
                protons,
                key=lambda idx: system.get_distance(donor, idx, mic=False),
            )
        ]
    return protons, acceptors


def group_acceptors_by_element(
    system: Atoms,
    acceptors: List[int],
) -> Dict[str, List[int]]:
    groups: Dict[str, List[int]] = {"O": [], "N": []}
    for idx in acceptors:
        symbol = system[idx].symbol
        if symbol in groups:
            groups[symbol].append(idx)
    return {key: value for key, value in groups.items() if value}


def _nearest_distance_lines(
    scope: ProtonTransferScope,
    protons: List[int],
    acceptors: List[int],
    idx_start_from: int,
    min_beta: float,
) -> tuple[List[str], str, List[str]]:
    donor_pl = scope.index_donor + 1 - idx_start_from
    protons_pl = _plumed_atom_list(protons, idx_start_from)
    acceptors_pl = _plumed_atom_list(acceptors, idx_start_from)
    beta = min_beta
    lines = [
        (
            f"pt_d_donor: DISTANCES GROUPA={donor_pl} GROUPB={protons_pl} "
            f"NOPBC MIN={{BETA={beta}}}"
        ),
        (
            f"pt_d_acc: DISTANCES GROUPA={acceptors_pl} GROUPB={protons_pl} "
            f"NOPBC MIN={{BETA={beta}}}"
        ),
        "pt_cv: COMBINE ARG=pt_d_acc.min,pt_d_donor.min COEFFICIENTS=1,-1 PERIODIC=NO",
    ]
    return lines, "pt_cv", ["pt_d_donor.min", "pt_d_acc.min"]


def _coordination_lines(
    scope: ProtonTransferScope,
    protons: List[int],
    acceptors: List[int],
    idx_start_from: int,
    cv_config: ProtonTransferCVConfig,
) -> tuple[List[str], str, List[str]]:
    if len(protons) != 1:
        raise ValueError("coordination mode requires a single transfer proton")
    proton_pl = protons[0] + 1 - idx_start_from
    donor_pl = scope.index_donor + 1 - idx_start_from
    acceptors_pl = _plumed_atom_list(acceptors, idx_start_from)
    r0 = cv_config.coordination_r0
    d0 = cv_config.coordination_d0
    lines = [
        (
            f"pt_cn_acc: COORDINATION GROUPA={proton_pl} GROUPB={acceptors_pl} "
            f"R_0={r0} D_0={d0} NN=6 MM=12 NOPBC"
        ),
        (
            f"pt_cn_don: COORDINATION GROUPA={proton_pl} GROUPB={donor_pl} "
            f"R_0={r0} D_0={d0} NN=6 MM=12 NOPBC"
        ),
        "pt_cv: COMBINE ARG=pt_cn_acc,pt_cn_don COEFFICIENTS=1,-1 PERIODIC=NO",
    ]
    return lines, "pt_cv", ["pt_cn_acc", "pt_cn_don"]


def _grouped_acceptor_lines(
    scope: ProtonTransferScope,
    protons: List[int],
    acceptors: List[int],
    idx_start_from: int,
    min_beta: float,
    system: Atoms,
    cv_config: ProtonTransferCVConfig,
) -> tuple[List[str], str, List[str]]:
    if len(protons) != 1:
        raise ValueError("grouped_acceptor mode requires a single transfer proton")
    donor_pl = scope.index_donor + 1 - idx_start_from
    proton_pl = protons[0] + 1 - idx_start_from
    beta = min_beta
    if cv_config.acceptor_grouping == "element":
        groups = group_acceptors_by_element(system, acceptors)
    else:
        groups = {"all": acceptors}
    if not groups:
        raise ValueError("grouped_acceptor mode has no acceptor groups")

    gap_names: List[str] = []
    print_extra: List[str] = []
    lines: List[str] = []
    for label, group_indices in groups.items():
        acceptors_pl = _plumed_atom_list(group_indices, idx_start_from)
        acc_label = label.lower()
        lines.append(f"pt_d_don_{acc_label}: DISTANCE ATOMS={donor_pl},{proton_pl} NOPBC")
        if len(group_indices) == 1:
            lines.append(
                f"pt_d_acc_{acc_label}: DISTANCE ATOMS={acceptors_pl},{proton_pl} NOPBC"
            )
            gap_name = f"pt_gap_{acc_label}"
            lines.append(
                f"{gap_name}: COMBINE ARG=pt_d_acc_{acc_label},pt_d_don_{acc_label} "
                "COEFFICIENTS=1,-1 PERIODIC=NO"
            )
        else:
            lines.append(
                (
                    f"pt_d_acc_{acc_label}: DISTANCES GROUPA={acceptors_pl} "
                    f"GROUPB={proton_pl} NOPBC MIN={{BETA={beta}}}"
                )
            )
            gap_name = f"pt_gap_{acc_label}"
            lines.append(
                f"{gap_name}: COMBINE ARG=pt_d_acc_{acc_label}.min,pt_d_don_{acc_label} "
                "COEFFICIENTS=1,-1 PERIODIC=NO"
            )
            print_extra.append(f"pt_d_acc_{acc_label}.min")
        gap_names.append(gap_name)
        print_extra.append(f"pt_d_don_{acc_label}")

    if len(gap_names) == 1:
        lines.append(f"pt_cv: COMBINE ARG={gap_names[0]} COEFFICIENTS=1 PERIODIC=NO")
    else:
        lines.append(f"pt_sort: SORT ARG={','.join(gap_names)}")
        lines.append("pt_cv: COMBINE ARG=pt_sort.1 COEFFICIENTS=1 PERIODIC=NO")
        print_extra.append("pt_sort.1")
    print_extra.extend(gap_names)
    return lines, "pt_cv", print_extra


def _coupled_dd_pt_lines(
    scope: ProtonTransferScope,
    protons: List[int],
    acceptors: List[int],
    idx_start_from: int,
    min_beta: float,
    cv_config: ProtonTransferCVConfig,
) -> tuple[List[str], str, List[str]]:
    base_lines, _, print_extra = _nearest_distance_lines(
        scope, protons, acceptors, idx_start_from, min_beta
    )
    donor_pl = scope.index_donor + 1 - idx_start_from
    acceptors_pl = _plumed_atom_list(acceptors, idx_start_from)
    beta = min_beta
    env_lines = [
        (
            f"pt_d_donacc: DISTANCES GROUPA={donor_pl} GROUPB={acceptors_pl} "
            f"NOPBC MIN={{BETA={beta}}}"
        ),
    ]
    lines = base_lines[:-1] + env_lines
    lines.append(
        "pt_core: COMBINE ARG=pt_d_acc.min,pt_d_donor.min COEFFICIENTS=1,-1 PERIODIC=NO"
    )
    lines.append(
        "pt_cv: COMBINE ARG=pt_core,pt_d_donacc.min,dd "
        f"COEFFICIENTS={cv_config.coupled_pt_weight},{cv_config.coupled_env_weight},"
        f"{cv_config.coupled_dd_weight} PERIODIC=NO"
    )
    print_extra.extend(["pt_core", "pt_d_donacc.min"])
    return lines, "pt_cv", print_extra


def generate_proton_transfer_geometry_walls(
    scope: ProtonTransferScope,
    protons: List[int],
    acceptors: List[int],
    idx_start_from: int,
    cv_config: ProtonTransferCVConfig,
    min_beta: float,
) -> List[str]:
    if not cv_config.geometry_walls:
        return []

    donor_pl = scope.index_donor + 1 - idx_start_from
    protons_pl = _plumed_atom_list(protons, idx_start_from)
    acceptors_pl = _plumed_atom_list(acceptors, idx_start_from)
    beta = min_beta
    kappa = 1000 * kcal / mol
    walls = [
        (
            f"pt_d_donor_wall: DISTANCES GROUPA={donor_pl} GROUPB={protons_pl} "
            f"NOPBC MIN={{BETA={beta}}}"
        ),
        (
            f"pt_uwall_h: UPPER_WALLS ARG=pt_d_donor_wall.min "
            f"AT={cv_config.max_donor_h_distance} KAPPA={kappa}"
        ),
        (
            f"pt_d_donacc_wall: DISTANCES GROUPA={donor_pl} GROUPB={acceptors_pl} "
            f"NOPBC MIN={{BETA={beta}}}"
        ),
        (
            f"pt_uwall_donacc: UPPER_WALLS ARG=pt_d_donacc_wall.min "
            f"AT={cv_config.max_donor_acceptor_distance} KAPPA={kappa}"
        ),
    ]
    return walls


def generate_proton_transfer_cv_lines(
    scope: ProtonTransferScope,
    idx_start_from: int,
    min_beta: float = 50.0,
    cv_config: Optional[ProtonTransferCVConfig] = None,
    system: Optional[Atoms] = None,
) -> tuple[List[str], str, List[str]]:
    if not scope.transfer_protons:
        raise ValueError("Proton transfer scope has no transfer protons")
    if not scope.acceptors:
        raise ValueError("Proton transfer scope has no acceptors")

    if cv_config is None:
        cv_config = ProtonTransferCVConfig()

    protons, acceptors = _resolve_scope_atoms(scope, cv_config, system)
    if not protons:
        raise ValueError("No transfer protons available for proton-transfer CV")
    if not acceptors:
        raise ValueError("No acceptors available for proton-transfer CV")

    mode = cv_config.mode
    if mode == "nearest_distance":
        lines, cv_name, print_extra = _nearest_distance_lines(
            scope, protons, acceptors, idx_start_from, min_beta
        )
    elif mode == "coordination":
        lines, cv_name, print_extra = _coordination_lines(
            scope, protons, acceptors, idx_start_from, cv_config
        )
    elif mode == "grouped_acceptor":
        if system is None:
            raise ValueError("grouped_acceptor mode requires the ASE system")
        lines, cv_name, print_extra = _grouped_acceptor_lines(
            scope,
            protons,
            acceptors,
            idx_start_from,
            min_beta,
            system,
            cv_config,
        )
    elif mode == "coupled_dd_pt":
        lines, cv_name, print_extra = _coupled_dd_pt_lines(
            scope, protons, acceptors, idx_start_from, min_beta, cv_config
        )
    else:
        raise ValueError(f"Unknown proton-transfer CV mode: {mode!r}")

    wall_lines = generate_proton_transfer_geometry_walls(
        scope, protons, acceptors, idx_start_from, cv_config, min_beta
    )
    return lines + wall_lines, cv_name, print_extra


def generate_opes_explore(
    cv_name: str = "pt_cv",
    pace: int = 20,
    barrier_kj_mol: float = 1.0,
    biasfactor: float = 15.0,
    state_wfile: str = "opes_state.data",
    state_wstride: int = 500,
    state_rfile: Optional[str] = None,
) -> tuple[List[str], str]:
    barrier = barrier_kj_mol * kJ / mol
    line = (
        f"opes: OPES_METAD_EXPLORE ARG={cv_name} PACE={pace} BARRIER={barrier} "
        f"BIASFACTOR={biasfactor} COMPRESSION_THRESHOLD=0.5 "
        f"STATE_WFILE={state_wfile} STATE_WSTRIDE={state_wstride}"
    )
    if state_rfile is not None:
        line += f" STATE_RFILE={state_rfile}"
    return [line], "opes.*"


def append_proton_transfer_to_plumed(
    plumed_config: List[str],
    scope: ProtonTransferScope,
    idx_start_from: int,
    dump_interval: int,
    bias: Optional[ProtonTransferBias] = None,
    cv_config: Optional[ProtonTransferCVConfig] = None,
    system: Optional[Atoms] = None,
) -> List[str]:
    if bias is None:
        bias = ProtonTransferBias()

    config = list(plumed_config)
    if bias.restart:
        if not config or config[0] != "RESTART":
            config.insert(0, "RESTART")

    cv_lines, cv_name, print_extra = generate_proton_transfer_cv_lines(
        scope,
        idx_start_from,
        min_beta=bias.min_beta,
        cv_config=cv_config,
        system=system,
    )
    state_rfile = bias.state_wfile if bias.restart else None
    opes_lines, opes_print = generate_opes_explore(
        cv_name=cv_name,
        pace=bias.pace,
        barrier_kj_mol=bias.barrier_kj_mol,
        biasfactor=bias.biasfactor,
        state_wfile=bias.state_wfile,
        state_wstride=bias.state_wstride,
        state_rfile=state_rfile,
    )

    print_idx = next(i for i, line in enumerate(config) if line.startswith("PRINT ARG="))
    for offset, line in enumerate(cv_lines + opes_lines):
        config.insert(print_idx + offset, line)

    print_idx = next(i for i, line in enumerate(config) if line.startswith("PRINT ARG="))
    print_line = config[print_idx]
    arg_part, stride_part = print_line.split(" STRIDE=", 1)
    old_args = arg_part.split("PRINT ARG=", 1)[1]
    extra = ",".join(dict.fromkeys(print_extra))
    merged_args = f"{old_args},{cv_name},{extra},{opes_print}" if extra else f"{old_args},{cv_name},{opes_print}"
    config[print_idx] = f"PRINT ARG={merged_args} STRIDE={stride_part}"
    return config


def _external_to_ase_index(index: int, idx_start_from: int) -> int:
    return int(index) - idx_start_from


def _selector_to_external_index(
    selector: AtomSelector,
    indices: Dict[str, int],
) -> int:
    if isinstance(selector, str):
        if selector not in indices:
            raise KeyError(f"Unknown atom selector {selector!r}")
        return int(indices[selector])
    return int(selector)


def _selectors_to_ase_indices(
    selectors: Optional[Union[AtomSelector, List[AtomSelector]]],
    indices: Dict[str, int],
    idx_start_from: int,
) -> Optional[List[int]]:
    values = _as_list(selectors)
    if values is None:
        return None
    return [
        _external_to_ase_index(_selector_to_external_index(value, indices), idx_start_from)
        for value in values
    ]


def _load_scope_file(scope_file: str) -> ProtonTransferScope:
    with open(scope_file, "r") as handle:
        return ProtonTransferScope.from_dict(json.load(handle))


def _save_scope_file(scope_file: str, scope: ProtonTransferScope) -> None:
    import os

    os.makedirs(os.path.dirname(os.path.abspath(scope_file)), exist_ok=True)
    with open(scope_file, "w") as handle:
        json.dump(scope.to_dict(), handle, indent=2)


def _covalent_radius_angstrom(symbol: str) -> float:
    try:
        from mendeleev import element
    except ImportError as exc:
        raise ImportError(
            "mendeleev is required for automatic proton-transfer detection"
        ) from exc
    radius_pm = element(symbol).covalent_radius
    if radius_pm is None:
        raise ValueError(f"No covalent radius available for element {symbol!r}")
    return radius_pm / 100.0


def _nearest_heavy_atom_index(system: Atoms, atom_idx: int) -> int:
    heavy_atoms = [i for i, atom in enumerate(system) if atom.number > 1]
    if not heavy_atoms:
        raise ValueError("No heavy atoms found in system")
    return min(heavy_atoms, key=lambda i: system.get_distance(atom_idx, i, mic=False))


def _heavy_bond_counts_from_mol(reference_mol, expected_n_atoms: int) -> Optional[List[int]]:
    if reference_mol is None:
        return None
    if reference_mol.GetNumAtoms() != expected_n_atoms:
        raise ValueError(
            f"Reference mol atom count ({reference_mol.GetNumAtoms()}) does not match "
            f"system atom count ({expected_n_atoms})"
        )
    counts = [0 for _ in range(expected_n_atoms)]
    for atom in reference_mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum() <= 1:
            continue
        counts[idx] = sum(
            1 for bond in atom.GetBonds() if bond.GetOtherAtom(atom).GetAtomicNum() > 1
        )
    return counts


def _is_topologically_allowed_acceptor(atomic_number: int, heavy_bond_count: int) -> bool:
    if atomic_number == 8:
        return heavy_bond_count < 2
    if atomic_number == 7:
        return heavy_bond_count < 3
    return False


class LocalOpesProtonTransferPlugin(ProtonTransferPlugin):
    """Local donor-H/acceptor proton-transfer CV with OPES exploration bias."""

    name = "local_opes"

    def _resolve_scope(
        self,
        generator: "PlumedConfigGenerator",
        config: ProtonTransferConfig,
    ) -> ProtonTransferScope:
        import os

        if config.scope_file is not None and os.path.exists(config.scope_file):
            return _load_scope_file(config.scope_file)

        indices = generator.get_indices()
        if config.donor is None:
            raise ValueError("proton-transfer donor was not configured")
        donor_external = _selector_to_external_index(config.donor, indices)
        donor = _external_to_ase_index(donor_external, generator.idx_start_from)
        if donor < 0 or donor >= len(generator.system):
            raise ValueError(f"donor index {donor_external} is out of range")

        explicit_protons = _selectors_to_ase_indices(
            config.proton, indices, generator.idx_start_from
        )
        explicit_acceptors = _selectors_to_ase_indices(
            config.acceptor, indices, generator.idx_start_from
        )
        if explicit_protons is not None and explicit_acceptors is not None:
            scope = ProtonTransferScope(
                index_donor=donor,
                transfer_protons=explicit_protons,
                acceptors=explicit_acceptors,
            )
            if config.scope_file is not None:
                _save_scope_file(config.scope_file, scope)
            return scope

        donor_symbol = generator.system[donor].symbol
        proton_threshold = (
            _covalent_radius_angstrom(donor_symbol) + _covalent_radius_angstrom("H")
        ) * config.proton_hbond_coeff
        transfer_protons = explicit_protons or [
            i
            for i, atom in enumerate(generator.system)
            if atom.number == 1
            and generator.system.get_distance(donor, i, mic=False) < proton_threshold
            and _nearest_heavy_atom_index(generator.system, i) == donor
        ]
        if not transfer_protons:
            raise ValueError(
                f"No transfer protons found within {proton_threshold:.3f} A of donor"
            )

        if explicit_acceptors is not None:
            acceptors = explicit_acceptors
        else:
            reference_mol = generator.reference_mol
            if reference_mol is None and config.topology_mol_file is not None:
                reference_mol = generator._load_reference_mol(config.topology_mol_file)
            heavy_bond_counts = _heavy_bond_counts_from_mol(
                reference_mol,
                len(generator.system),
            )
            candidates = [
                i
                for i, atom in enumerate(generator.system)
                if atom.number in (7, 8)
                and i != donor
                and (
                    heavy_bond_counts is None
                    or _is_topologically_allowed_acceptor(atom.number, heavy_bond_counts[i])
                )
            ]
            acceptors = [
                i
                for i in candidates
                if generator.system.get_distance(donor, i, mic=False) < config.acceptor_radius
            ]
        if not acceptors:
            raise ValueError(
                f"No acceptors found within {config.acceptor_radius:.3f} A of donor"
            )

        scope = ProtonTransferScope(
            index_donor=donor,
            transfer_protons=transfer_protons,
            acceptors=acceptors,
        )
        if config.scope_file is not None:
            _save_scope_file(config.scope_file, scope)
        return scope

    def append_to_plumed(
        self,
        generator: "PlumedConfigGenerator",
        plumed_config: List[str],
        config: ProtonTransferConfig,
        dump_interval: int,
        integrate_config: dict,
    ) -> List[str]:
        scope = self._resolve_scope(generator, config)
        state_wstride = (
            config.opes_state_wstride
            if config.opes_state_wstride is not None
            else default_opes_state_wstride(integrate_config.get("n_step"), dump_interval)
        )
        bias = ProtonTransferBias(
            pace=config.opes_pace,
            barrier_kj_mol=config.opes_barrier,
            biasfactor=config.opes_biasfactor,
            state_wfile=config.state_file,
            state_wstride=state_wstride,
            state_rfile=config.state_file if config.restart else None,
            restart=config.restart,
        )
        cv_config = ProtonTransferCVConfig(
            mode=config.flavor or "nearest_distance",
            acceptor_grouping=config.acceptor_grouping,
            coordination_r0=config.coordination_r0,
            coordination_d0=config.coordination_d0,
            coupled_pt_weight=config.coupled_pt_weight,
            coupled_env_weight=config.coupled_env_weight,
            coupled_dd_weight=config.coupled_dd_weight,
            geometry_walls=config.geometry_walls,
            max_donor_h_distance=config.max_donor_h_distance,
            max_donor_acceptor_distance=config.max_donor_acceptor_distance,
            single_proton=config.single_proton,
        )
        return append_proton_transfer_to_plumed(
            plumed_config,
            scope,
            generator.idx_start_from,
            dump_interval,
            bias=bias,
            cv_config=cv_config,
            system=generator.system,
        )


register_proton_transfer_plugin(LocalOpesProtonTransferPlugin())
