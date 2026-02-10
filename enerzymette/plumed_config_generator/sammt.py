from typing import List, Optional
from ase import Atoms
from ase.units import kJ, mol, fs, kcal

def get_sammt_basics(
    system: Atoms,
    idx_start_from: int,
    dump_interval: int,
    reference_pdb_file: Optional[str]=None,
    substrate: Optional[str]=None,
    nucleophile: Optional[str]=None,
    index_sulphur: Optional[int]=None,
    index_methyl_carbon: Optional[int]=None,
    index_nucleophile: Optional[int]=None,
    **kwargs
) -> List[str]:
    plumed_config = []
    if reference_pdb_file is not None:
        idx_start_from = 1
        with open(reference_pdb_file, "r") as f:
            atom_count = 0
            for line in f.readlines():
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    atom_count += 1
                    resname = line[17:20]
                    atomname = line[11:16].strip()
                    if resname == "SAM":
                        if atomname == "SD":
                            index_sulphur = atom_count
                        elif atomname == "CE":
                            index_methyl_carbon = atom_count
                    elif resname == substrate:
                        if atomname == nucleophile:
                            index_nucleophile = atom_count
    if index_sulphur is None or index_methyl_carbon is None or index_nucleophile is None:
        raise ValueError("Index of sulphur, methyl carbon, and nucleophile must be provided")
    plumed_config.append(f"UNITS LENGTH=A TIME={0.001 / fs} ENERGY={1 / kJ * mol}")
    plumed_config.append(f"d0: DISTANCE ATOMS={index_sulphur + 1 - idx_start_from},{index_methyl_carbon + 1 - idx_start_from} NOPBC")
    plumed_config.append(f"d1: DISTANCE ATOMS={index_methyl_carbon + 1 - idx_start_from},{index_nucleophile + 1 - idx_start_from} NOPBC")
    plumed_config.append(f"dsort: SORT ARG=d0,d1")
    plumed_config.append(f"uwall: UPPER_WALLS ARG=dsort.1 AT=3.0 KAPPA={1000 * kcal / mol}")
    plumed_config.append(f"dd: COMBINE ARG=d1,d0 COEFFICIENTS=1,-1 PERIODIC=NO")
    current_d0 = system.get_distance(index_sulphur - idx_start_from, index_methyl_carbon - idx_start_from)
    current_d1 = system.get_distance(index_methyl_carbon - idx_start_from, index_nucleophile - idx_start_from)
    current_dd = current_d1 - current_d0
    return plumed_config, current_dd


def get_sammt_config(
    system: Atoms,
    integrate_config: dict,
    idx_start_from: int,
    dump_interval: int,
    upper_bound: float,
    lower_bound: float,
    reference_pdb_file: Optional[str]=None,
    substrate: Optional[str]=None,
    nucleophile: Optional[str]=None,
    index_sulphur: Optional[int]=None,
    index_methyl_carbon: Optional[int]=None,
    index_nucleophile: Optional[int]=None,
    **kwargs
) -> List[str]:
    plumed_config, current_dd = get_sammt_basics(system, idx_start_from, dump_interval, reference_pdb_file, substrate, nucleophile, index_sulphur, index_methyl_carbon, index_nucleophile)
    
    interval = upper_bound - lower_bound
    n_step = integrate_config.get("n_step")
    moving_restraint_kappa = 1000 * kcal / mol
    moving_restraint_component = [f"mr: MOVINGRESTRAINT ARG=dd STEP0=0 AT0={current_dd} KAPPA0={moving_restraint_kappa}"]
    if current_dd <= lower_bound:
        step1 = n_step // 2
        stage = len(moving_restraint_component)
        moving_restraint_component.append(f"STEP{stage}={step1} AT{stage}={upper_bound}")
        stop_dd = lower_bound
    elif current_dd <= upper_bound:
        step1 = int(n_step // 2 * (upper_bound - current_dd) / interval)
        if step1 > 0:
            stage = len(moving_restraint_component)
            moving_restraint_component.append(f"STEP{stage}={step1} AT{stage}={upper_bound}")
            stop_dd = current_dd
        else:
            stop_dd = upper_bound
    else:
        step1 = 0
        stop_dd = upper_bound
    step2 = n_step // 2 + step1
    stage = len(moving_restraint_component)
    moving_restraint_component.append(f"STEP{stage}={step2} AT{stage}={lower_bound}")
    if n_step - step2 > 0:
        stage = len(moving_restraint_component)
        moving_restraint_component.append(f"STEP{stage}={n_step} AT{stage}={stop_dd}")
    plumed_config.append(" ".join(moving_restraint_component))
    plumed_config.append(f"PRINT ARG=d0,d1,dsort.1,dd,mr.* STRIDE={dump_interval}")
    plumed_config.append(f"FLUSH STRIDE={dump_interval}")
    return plumed_config

def get_naive_sammt_config(
    system: Atoms,
    integrate_config: dict,
    idx_start_from: int,
    dump_interval: int,
    upper_bound: float,
    lower_bound: float,
    warmup_steps: int,
    reference_pdb_file: Optional[str]=None,
    substrate: Optional[str]=None,
    nucleophile: Optional[str]=None,
    index_sulphur: Optional[int]=None,
    index_methyl_carbon: Optional[int]=None,
    index_nucleophile: Optional[int]=None,
    **kwargs
) -> List[str]:
    plumed_config, current_dd = get_sammt_basics(system, idx_start_from, dump_interval, reference_pdb_file, substrate, nucleophile, index_sulphur, index_methyl_carbon, index_nucleophile)
    n_step = integrate_config.get("n_step")
    moving_restraint_kappa = 1000 * kcal / mol
    moving_restraint_component = f"mr: MOVINGRESTRAINT ARG=dd STEP0=0 AT0={current_dd} KAPPA0={moving_restraint_kappa} STEP1={warmup_steps} AT1={upper_bound} STEP2={n_step} AT2={lower_bound}"
    plumed_config.append(moving_restraint_component)
    plumed_config.append(f"PRINT ARG=d0,d1,dsort.1,dd,mr.* STRIDE={dump_interval}")
    plumed_config.append(f"FLUSH STRIDE={dump_interval}")
    return plumed_config
