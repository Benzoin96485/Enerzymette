from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union
import warnings

from ase import Atoms
from ase.units import fs, kJ, kcal, mol

from .proton_transfer import (
    ProtonTransferConfig,
    append_optional_proton_transfer,
    build_proton_transfer_config,
)


@dataclass
class ReactionCoordinate:
    preamble: List[str]
    cv_name: str
    lower_bound: float
    upper_bound: float
    initial_value: float
    dump_interval: int
    kappa: float = field(default_factory=lambda: 1000 * kcal / mol)
    print_args: Optional[str] = None


class PlumedConfigGenerator(ABC):
    """Abstract base class for dynamic PLUMED config generators.

    Subclasses provide chemistry-specific atom indexing and the main reaction
    coordinate. This base class owns generic PLUMED workflows and optional
    proton-transfer plugin insertion.
    """

    default_cv_name: str = "rc"
    default_print_args: Optional[str] = None

    def __init__(
        self,
        system: Atoms,
        *,
        integrate_config: Optional[dict] = None,
        idx_start_from: int = 1,
        preamble: Optional[List[str]] = None,
        reference_pdb: Optional[str] = None,
        reference_pdb_file: Optional[str] = None,
        reference_mol=None,
        reference_mol_file: Optional[str] = None,
        topology_mol_file: Optional[str] = None,
        proton_transfer: Optional[Union[bool, dict, ProtonTransferConfig]] = None,
        proton_transfer_plugin: Optional[str] = None,
        **kwargs,
    ) -> None:
        self.system = system
        self.integrate_config = integrate_config or {}
        self.idx_start_from = idx_start_from
        self.preamble = preamble or [
            f"UNITS LENGTH=A TIME={0.001 / fs} ENERGY={1 / kJ * mol}"
        ]
        self.reference_pdb = reference_pdb_file or reference_pdb
        self.reference_mol = reference_mol
        self.reference_mol_file = reference_mol_file or topology_mol_file
        if self.reference_mol is None and self.reference_mol_file is not None:
            self.reference_mol = self._load_reference_mol(self.reference_mol_file)
        self.proton_transfer_config = build_proton_transfer_config(
            proton_transfer,
            proton_transfer_plugin=proton_transfer_plugin,
            topology_mol_file=topology_mol_file,
        )

    @staticmethod
    def _load_reference_mol(path: str):
        try:
            from rdkit import Chem
        except ImportError:
            warnings.warn(
                "RDKit is not available; reference_mol was not loaded.",
                RuntimeWarning,
            )
            return None
        mol = Chem.MolFromMolFile(path, removeHs=False, sanitize=False)
        if mol is None:
            warnings.warn(f"Could not load reference mol from {path}", RuntimeWarning)
        return mol

    @abstractmethod
    def get_indices(self) -> Dict[str, int]:
        """Return descriptive atom names mapped to atom indices."""

    @abstractmethod
    def define_main_rc(self) -> Tuple[str, str]:
        """Return the main CV name and its PLUMED definition lines."""

    @abstractmethod
    def calc_main_rc(self) -> float:
        """Calculate the current main CV value from ``self.system``."""

    def build_reaction_coordinate(
        self,
        *,
        lower_bound: float,
        upper_bound: float,
        dump_interval: int,
        kappa: Optional[float] = None,
        print_args: Optional[str] = None,
        **kwargs,
    ) -> ReactionCoordinate:
        cv_name, definition = self.define_main_rc()
        lines = list(self.preamble)
        lines.extend(line for line in definition.splitlines() if line.strip())
        return ReactionCoordinate(
            preamble=lines,
            cv_name=cv_name,
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            initial_value=self.calc_main_rc(),
            dump_interval=dump_interval,
            kappa=1000 * kcal / mol if kappa is None else kappa,
            print_args=print_args if print_args is not None else self.default_print_args,
        )

    def _append_optional_proton_transfer(
        self,
        plumed_config: List[str],
        *,
        dump_interval: int,
        integrate_config: Optional[dict] = None,
        proton_transfer: Optional[Union[bool, dict, ProtonTransferConfig]] = None,
        proton_transfer_plugin: Optional[str] = None,
        **kwargs,
    ) -> List[str]:
        return append_optional_proton_transfer(
            self,
            plumed_config,
            dump_interval=dump_interval,
            integrate_config=integrate_config,
            proton_transfer=proton_transfer,
            proton_transfer_plugin=proton_transfer_plugin,
        )

    def standard_steered_md(
        self,
        *,
        integrate_config: Optional[dict] = None,
        lower_bound: float,
        upper_bound: float,
        dump_interval: int,
        **kwargs,
    ) -> List[str]:
        rc = self.build_reaction_coordinate(
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            dump_interval=dump_interval,
            **kwargs,
        )
        plumed_config = generate_steered_md(rc, integrate_config or self.integrate_config)
        return self._append_optional_proton_transfer(
            plumed_config,
            dump_interval=dump_interval,
            integrate_config=integrate_config,
            **kwargs,
        )

    def naive_steered_md(
        self,
        *,
        integrate_config: Optional[dict] = None,
        lower_bound: float,
        upper_bound: float,
        dump_interval: int,
        warmup_steps: int,
        **kwargs,
    ) -> List[str]:
        rc = self.build_reaction_coordinate(
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            dump_interval=dump_interval,
            **kwargs,
        )
        n_step = (integrate_config or self.integrate_config).get("n_step")
        if n_step is None:
            raise ValueError("integrate_config.n_step is required for naive steered MD")
        lower_dist = abs(rc.initial_value - lower_bound)
        upper_dist = abs(rc.initial_value - upper_bound)
        first, second = (
            (lower_bound, upper_bound) if lower_dist <= upper_dist else (upper_bound, lower_bound)
        )
        plumed_config = list(rc.preamble)
        plumed_config.append(
            f"mr: MOVINGRESTRAINT ARG={rc.cv_name} STEP0=0 AT0={rc.initial_value} "
            f"KAPPA0={rc.kappa} STEP1={warmup_steps} AT1={first} STEP2={n_step} AT2={second}"
        )
        print_args = rc.print_args if rc.print_args is not None else f"{rc.cv_name},mr.*"
        plumed_config.append(f"PRINT ARG={print_args} STRIDE={dump_interval}")
        plumed_config.append(f"FLUSH STRIDE={dump_interval}")
        return self._append_optional_proton_transfer(
            plumed_config,
            dump_interval=dump_interval,
            integrate_config=integrate_config,
            **kwargs,
        )

    def scan(
        self,
        *,
        target_value: float,
        lower_bound: float,
        upper_bound: float,
        dump_interval: int,
        **kwargs,
    ) -> List[str]:
        rc = self.build_reaction_coordinate(
            lower_bound=lower_bound,
            upper_bound=upper_bound,
            dump_interval=dump_interval,
            **kwargs,
        )
        return generate_scan_restraint(rc, target_value)




def generate_steered_md(rc: ReactionCoordinate, integrate_config: dict) -> List[str]:
    plumed_config = list(rc.preamble)

    interval = rc.upper_bound - rc.lower_bound
    n_step = integrate_config.get("n_step")
    moving_restraint_kappa = rc.kappa
    moving_restraint_component = [
        f"mr: MOVINGRESTRAINT ARG={rc.cv_name} STEP0=0 AT0={rc.initial_value} KAPPA0={moving_restraint_kappa}"
    ]
    if rc.initial_value <= rc.lower_bound:
        step1 = n_step // 2
        stage = len(moving_restraint_component)
        moving_restraint_component.append(f"STEP{stage}={step1} AT{stage}={rc.upper_bound}")
        stop_dd = rc.lower_bound
    elif rc.initial_value <= rc.upper_bound:
        step1 = int(n_step // 2 * (rc.upper_bound - rc.initial_value) / interval)
        if step1 > 0:
            stage = len(moving_restraint_component)
            moving_restraint_component.append(f"STEP{stage}={step1} AT{stage}={rc.upper_bound}")
            stop_dd = rc.initial_value
        else:
            stop_dd = rc.upper_bound
    else:
        step1 = 0
        stop_dd = rc.upper_bound
    step2 = n_step // 2 + step1
    stage = len(moving_restraint_component)
    moving_restraint_component.append(f"STEP{stage}={step2} AT{stage}={rc.lower_bound}")
    if n_step - step2 > 0:
        stage = len(moving_restraint_component)
        moving_restraint_component.append(f"STEP{stage}={n_step} AT{stage}={stop_dd}")
    plumed_config.append(" ".join(moving_restraint_component))

    print_args = rc.print_args if rc.print_args is not None else f"{rc.cv_name},mr.*"
    plumed_config.append(f"PRINT ARG={print_args} STRIDE={rc.dump_interval}")
    plumed_config.append(f"FLUSH STRIDE={rc.dump_interval}")
    return plumed_config


def generate_scan_restraint(
    rc: ReactionCoordinate,
    target_value: float,
    kappa: Optional[float] = None,
) -> List[str]:
    plumed_config = list(rc.preamble)
    restraint_kappa = rc.kappa if kappa is None else kappa
    plumed_config.append(
        f"r: RESTRAINT ARG={rc.cv_name} AT={target_value} KAPPA={restraint_kappa}"
    )

    if rc.print_args is not None:
        print_args = rc.print_args.replace("mr.*", "r.*")
    else:
        print_args = f"{rc.cv_name},r.*"
    plumed_config.append(f"PRINT ARG={print_args} STRIDE={rc.dump_interval}")
    plumed_config.append(f"FLUSH STRIDE={rc.dump_interval}")
    return plumed_config
