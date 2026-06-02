from dataclasses import dataclass, field
from typing import List, Optional

from ase.units import kcal, mol


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
