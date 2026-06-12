from dataclasses import dataclass, field
from typing import List, Optional

from ase.units import kJ, kcal, mol


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
    pt_restart: bool = False
    min_beta: float = 50.0


def _plumed_atom_list(indices: List[int], idx_start_from: int) -> str:
    return ",".join(str(i + 1 - idx_start_from) for i in indices)


def generate_proton_transfer_cv_lines(
    scope: ProtonTransferScope,
    idx_start_from: int,
    min_beta: float = 50.0,
) -> tuple[List[str], str]:
    if not scope.transfer_protons:
        raise ValueError("Proton transfer scope has no transfer protons")
    if not scope.acceptors:
        raise ValueError("Proton transfer scope has no acceptors")

    donor_pl = scope.index_donor + 1 - idx_start_from
    protons_pl = _plumed_atom_list(scope.transfer_protons, idx_start_from)
    acceptors_pl = _plumed_atom_list(scope.acceptors, idx_start_from)
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
    return lines, "pt_cv"


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
) -> List[str]:
    if bias is None:
        bias = ProtonTransferBias()

    config = list(plumed_config)
    if bias.pt_restart:
        if not config or config[0] != "RESTART":
            config.insert(0, "RESTART")

    cv_lines, cv_name = generate_proton_transfer_cv_lines(
        scope, idx_start_from, min_beta=bias.min_beta
    )
    state_rfile = bias.state_wfile if bias.pt_restart else None
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
    config[print_idx] = (
        f"PRINT ARG={old_args},{cv_name},{opes_print} STRIDE={stride_part}"
    )
    return config


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
