"""Shared scan workflow helpers used by scantoolkit and altoolkit."""

from __future__ import annotations

import copy
import json
import os
import shutil
import subprocess
import sys
import traceback
from glob import glob
from typing import Any, Callable, Dict, List, Optional, Tuple

import ase.io
import yaml
from ase.units import kcal, mol

from ..logger import logger
from ..mep_util import analyze_scan_path, find_new_name
from ..plumed_config_generator import (
    get_config_generator_name,
    get_scan_method_name,
    resolve_scan_endpoints,
)

WriteConfigFn = Callable[..., None]
SimulateCmdFn = Callable[[str, str, bool], List[str]]


def build_enerzyme_simulate_cmd(
    config_path: str,
    output_path: str,
    model_path: str,
    *,
    calculator_patch: Optional[str] = None,
    plumed_patch: Optional[str] = None,
    model_config_path: Optional[str] = None,
    model_config_arg: Optional[List[str]] = None,
) -> List[str]:
    cmd = [
        "enerzyme", "simulate",
        "-c", config_path,
        "-o", output_path,
        "-m", model_path,
    ]
    if calculator_patch is not None:
        cmd.extend(["-cp", calculator_patch])
    if model_config_path is not None:
        cmd.extend(["-mc", model_config_path])
    elif model_config_arg:
        cmd.extend(model_config_arg)
    if plumed_patch is not None:
        cmd.extend(["-pp", plumed_patch])
    return cmd


def apply_plumed_scan_sampling(
    config: Dict[str, Any],
    *,
    structure_path: str,
    plumed_patch_key: str,
    plumed_cv_config: Dict[str, Any],
    n_steps: int,
    target_value: Optional[float] = None,
    target_structure_path: Optional[str] = None,
) -> None:
    config["Simulation"].pop("integrate", None)
    if "optimize" not in config["Simulation"]:
        config["Simulation"]["optimize"] = {"optimizer": "LBFGS"}
    idx_start_from = config["Simulation"]["idx_start_from"]
    initial_structure = ase.io.read(structure_path, index=-1)
    x0, x1, num, rc = resolve_scan_endpoints(
        initial_structure,
        idx_start_from,
        plumed_patch_key,
        plumed_cv_config,
        n_steps,
        target_value=target_value,
        target_structure_path=target_structure_path,
    )
    config["Simulation"]["plumed_config_generator"] = {
        "name": get_config_generator_name(plumed_patch_key),
        "method": get_scan_method_name(plumed_patch_key),
    }
    config["Simulation"]["sampling"] = {
        "cv": "plumed",
        "params": {
            "x0": float(x0),
            "x1": float(x1),
            "num": num,
            "plumed_config": dict(plumed_cv_config),
        },
    }
    logger.info(
        f"PLUMED scan on CV {rc.cv_name} from {x0} to {x1} "
        f"(bounds [{rc.lower_bound}, {rc.upper_bound}]) with {num} steps"
    )


def apply_bond_scan_sampling(
    config: Dict[str, Any],
    *,
    structure_path: str,
    constraint_scan: Dict[str, Any],
    n_steps: int,
    idx_start_from: int,
    target_value: Optional[float] = None,
    target_structure_path: Optional[str] = None,
) -> None:
    for constraint_type, constraint_params in constraint_scan.items():
        if constraint_type != "bond":
            continue
        index0 = constraint_params["i0"]
        index1 = constraint_params["i1"]
        ase_i0 = index0 - idx_start_from
        ase_i1 = index1 - idx_start_from
        initial_structure = ase.io.read(structure_path, index=-1)
        initial_value = initial_structure.get_distance(ase_i0, ase_i1)
        if target_value is None:
            if target_structure_path is None:
                from mendeleev import element

                element0 = element(initial_structure.symbols[ase_i0])
                element1 = element(initial_structure.symbols[ase_i1])
                target_value = (
                    element0.covalent_radius_pyykko + element1.covalent_radius_pyykko
                ) / 100
                logger.info(
                    f"Target value is set to {target_value} Angstrom based on "
                    f"single-bond Pyykko covalent radii for bond between atoms "
                    f"{index0} and {index1}"
                )
            else:
                target_structure = ase.io.read(target_structure_path, index=-1)
                target_value = target_structure.get_distance(ase_i0, ase_i1)
                logger.info(
                    f"Target value is set to {target_value} Angstrom based on "
                    f"distance between atoms {index0} and {index1} in reference "
                    f"target structure {target_structure_path}"
                )
        config["Simulation"]["sampling"] = {
            "cv": "distance",
            "params": {
                "x0": float(initial_value),
                "x1": float(target_value),
                "num": n_steps,
                "i0": index0,
                "i1": index1,
            },
        }
        logger.info(
            f"Scanning distance between atoms {index0} and {index1} from "
            f"{initial_value} to {target_value} with {n_steps} steps"
        )


def write_standalone_scan_config(
    config_path: str,
    *,
    task: str,
    initial_structure_path: str,
    charge: int,
    multiplicity: int,
    constraint_freeze_xyz: List[int],
    idx_start_from: int,
    constraint_scan: Optional[Dict[str, Any]] = None,
    plumed_patch_key: Optional[str] = None,
    plumed_cv_config: Optional[Dict[str, Any]] = None,
    n_steps: int = 25,
    target_value: Optional[float] = None,
    target_structure_path: Optional[str] = None,
    traj_file: Optional[str] = None,
) -> None:
    optimize: Dict[str, Any] = {
        "optimizer": "LBFGS",
    }
    if traj_file is not None:
        optimize["traj_file"] = traj_file
    base_config: Dict[str, Any] = {
        "Simulation": {
            "environment": "ase",
            "dtype": "float64",
            "cuda": True,
            "task": task,
            "idx_start_from": idx_start_from,
            "neighbor_list": "full",
            "constraint": {
                "fix_atom": {
                    "indices": constraint_freeze_xyz,
                }
            },
            "optimize": optimize,
        },
        "System": {
            "structure_file": initial_structure_path,
            "charge": charge,
            "multiplicity": multiplicity,
        },
    }
    if task == "plumed_scan":
        if plumed_patch_key is None:
            raise ValueError("plumed_scan task requires plumed_patch_key")
        apply_plumed_scan_sampling(
            base_config,
            structure_path=initial_structure_path,
            plumed_patch_key=plumed_patch_key,
            plumed_cv_config=plumed_cv_config or {},
            n_steps=n_steps,
            target_value=target_value,
            target_structure_path=target_structure_path,
        )
    elif task == "scan":
        apply_bond_scan_sampling(
            base_config,
            structure_path=initial_structure_path,
            constraint_scan=constraint_scan or {},
            n_steps=n_steps,
            idx_start_from=idx_start_from,
            target_value=target_value,
            target_structure_path=target_structure_path,
        )

    with open(config_path, "w") as handle:
        yaml.dump(base_config, handle)
    logger.info(f"Config (task: {task}) written to {config_path}")


def write_base_scan_task_config(
    config_path: str,
    *,
    task: str,
    structure_path: str,
    base_config: Dict[str, Any],
    plumed_patch_key: Optional[str] = None,
    n_steps: int = 25,
    scan_target_value: Optional[float] = None,
    scan_target_structure_path: Optional[str] = None,
    traj_file: Optional[str] = None,
    yaml_flow_style: bool = False,
) -> None:
    config = copy.deepcopy(base_config)
    config["System"]["structure_file"] = structure_path
    config["Simulation"]["task"] = task
    config["Simulation"].pop("uncertainty_calculator", None)

    if task == "opt":
        config["Simulation"].pop("sampling", None)
        config["Simulation"].pop("integrate", None)
        if "optimize" not in config["Simulation"]:
            config["Simulation"]["optimize"] = {"optimizer": "LBFGS"}
        if traj_file is not None:
            config["Simulation"]["optimize"]["traj_file"] = traj_file
    elif task == "plumed_scan":
        if plumed_patch_key is None:
            raise ValueError("plumed_scan task requires plumed_patch_key")
        plumed_cv_config = config["Simulation"]["sampling"]["params"]["plumed_config"]
        apply_plumed_scan_sampling(
            config,
            structure_path=structure_path,
            plumed_patch_key=plumed_patch_key,
            plumed_cv_config=plumed_cv_config,
            n_steps=n_steps,
            target_value=scan_target_value,
            target_structure_path=scan_target_structure_path,
        )

    with open(config_path, "w") as handle:
        yaml.dump(config, handle, default_flow_style=yaml_flow_style)


def copy_reaction_local_minima(
    local_minima_path: str,
    reactant_name: str,
    product_name: str,
    elementary_reaction_path: str,
    *,
    warn_on_missing: bool = True,
) -> None:
    reactant_src = os.path.join(elementary_reaction_path, "reactant.xyz")
    product_src = os.path.join(elementary_reaction_path, "product.xyz")
    if os.path.exists(reactant_src):
        shutil.copy(
            reactant_src,
            os.path.join(local_minima_path, "reactant", f"{reactant_name}.xyz"),
        )
    elif warn_on_missing:
        logger.warning(
            f"Reactant file does not exist for finished reaction "
            f"{reactant_name}-{product_name}"
        )
    if os.path.exists(product_src):
        shutil.copy(
            product_src,
            os.path.join(local_minima_path, "product", f"{product_name}.xyz"),
        )
    elif warn_on_missing:
        logger.warning(
            f"Product file does not exist for finished reaction "
            f"{reactant_name}-{product_name}"
        )


def find_lowest_local_minima(local_minima_path: str) -> Tuple[str, str, float, float]:
    reactant_paths = glob(os.path.join(local_minima_path, "reactant", "*.xyz"))
    product_paths = glob(os.path.join(local_minima_path, "product", "*.xyz"))
    reactant_energies = [
        ase.io.read(reactant_path, format="extxyz", index=-1).get_potential_energy()
        for reactant_path in reactant_paths
    ]
    product_energies = [
        ase.io.read(product_path, format="extxyz", index=-1).get_potential_energy()
        for product_path in product_paths
    ]
    lowest_reactant_energy = min(reactant_energies)
    lowest_product_energy = min(product_energies)
    lowest_reactant_name = os.path.basename(
        reactant_paths[reactant_energies.index(lowest_reactant_energy)]
    ).split(".")[0]
    lowest_product_name = os.path.basename(
        product_paths[product_energies.index(lowest_product_energy)]
    ).split(".")[0]
    return (
        lowest_reactant_name,
        lowest_product_name,
        lowest_reactant_energy,
        lowest_product_energy,
    )


def _run_simulate(cmd: List[str], *, redirect_stdio: bool = False) -> None:
    kwargs = {}
    if redirect_stdio:
        kwargs["stdout"] = sys.stdout
        kwargs["stderr"] = sys.stderr
    subprocess.Popen(cmd, **kwargs).wait()


def run_elementary_reaction_scan(
    reactant_path: str,
    elementary_reaction_path: str,
    *,
    write_config: WriteConfigFn,
    build_simulate_cmd: SimulateCmdFn,
    scan_task: str = "plumed_scan",
    target_value: Optional[float] = None,
    target_structure_path: Optional[str] = None,
    redirect_stdio: bool = False,
) -> Dict[str, Any]:
    if os.path.exists(elementary_reaction_path):
        logger.warning(
            f"Output path {elementary_reaction_path} already exists, it could be overwritten"
        )
    else:
        os.makedirs(elementary_reaction_path)

    logger.info(f"Reading reactant from {reactant_path}")
    reactant_atoms = ase.io.read(reactant_path, index=-1)
    init_reactant_path = os.path.join(elementary_reaction_path, "init_reactant.xyz")
    ase.io.write(init_reactant_path, reactant_atoms, format="extxyz")
    logger.info(f"Initial reactant written to {init_reactant_path}")

    reactant_opt_config_path = os.path.join(elementary_reaction_path, "reactant_opt.yaml")
    write_config(
        "opt",
        init_reactant_path,
        reactant_opt_config_path,
        traj_file="reactant_traj-opt.xyz",
    )
    _run_simulate(
        build_simulate_cmd(reactant_opt_config_path, elementary_reaction_path, False),
        redirect_stdio=redirect_stdio,
    )
    opt_path = os.path.join(elementary_reaction_path, "optim.xyz")
    optimized_reactant_path = os.path.join(elementary_reaction_path, "reactant.xyz")
    os.rename(opt_path, optimized_reactant_path)
    logger.info(f"Reactant optimized and written to {optimized_reactant_path}")

    scan_config_path = os.path.join(elementary_reaction_path, "scan.yaml")
    write_config(
        scan_task,
        optimized_reactant_path,
        scan_config_path,
        target_value=target_value,
        target_structure_path=target_structure_path,
    )
    _run_simulate(
        build_simulate_cmd(
            scan_config_path,
            elementary_reaction_path,
            scan_task == "plumed_scan",
        ),
        redirect_stdio=redirect_stdio,
    )
    logger.info(f"Scan finished for elementary reaction {elementary_reaction_path}")

    scan_atoms = ase.io.read(
        os.path.join(elementary_reaction_path, "scan_optim.xyz"),
        format="extxyz",
        index=":",
    )
    energies = [atoms.get_potential_energy() for atoms in scan_atoms]
    path = analyze_scan_path(energies)

    scan_product_path = os.path.join(elementary_reaction_path, "init_product.xyz")
    ase.io.write(scan_product_path, scan_atoms[path.product_index], format="extxyz")
    logger.info(
        f"Scanned product (image {path.product_index}) written to {scan_product_path}"
    )

    product_opt_config_path = os.path.join(elementary_reaction_path, "product_opt.yaml")
    write_config(
        "opt",
        scan_product_path,
        product_opt_config_path,
        traj_file="product_traj-opt.xyz",
    )
    _run_simulate(
        build_simulate_cmd(product_opt_config_path, elementary_reaction_path, False),
        redirect_stdio=redirect_stdio,
    )
    product_path = os.path.join(elementary_reaction_path, "product.xyz")
    os.rename(opt_path, product_path)
    logger.info(f"Product optimized and written to {product_path}")

    return {
        "intermediate_indices": path.intermediate_indices,
        "mep_path_info": {
            "atoms": scan_atoms,
            "energies": energies,
            "n_images": len(scan_atoms),
            "n_atoms": len(scan_atoms[0]),
        },
        "ci_index": path.ci_index,
        "terminate_chain": path.terminate_chain,
        "chain_reactant_index": path.chain_reactant_index,
    }


def write_rate_determining_ts_results(
    output_path: str,
    reaction_name: str,
    mep_path_info: Dict[str, Any],
    ci_index: int,
    local_minima_path: str,
) -> None:
    (
        lowest_reactant_name,
        lowest_product_name,
        lowest_reactant_energy,
        lowest_product_energy,
    ) = find_lowest_local_minima(local_minima_path)
    ts_energy = mep_path_info["energies"][ci_index]
    energy_span = (ts_energy - lowest_reactant_energy) / (kcal / mol)
    energy_change = (lowest_product_energy - lowest_reactant_energy) / (kcal / mol)
    logger.info(f"Reaction energy span: {energy_span:.2f} kcal/mol")
    logger.info(f"Reaction energy change: {energy_change:.2f} kcal/mol")
    results = {
        "energy span": energy_span,
        "energy change": energy_change,
        "lowest energy reactant": lowest_reactant_name,
        "lowest energy product": lowest_product_name,
    }
    ts_dir = os.path.join(output_path, "rate_determining_ts")
    os.makedirs(ts_dir, exist_ok=True)
    with open(os.path.join(ts_dir, "results.json"), "w") as handle:
        json.dump(results, handle, indent=4)
    ase.io.write(
        os.path.join(ts_dir, f"{reaction_name}.xyz"),
        mep_path_info["atoms"][ci_index],
        format="extxyz",
    )


def run_scan_chain(
    output_path: str,
    initial_reactant_path: str,
    *,
    run_elementary_reaction: Callable[..., Dict[str, Any]],
    local_minima_path: Optional[str] = None,
    reactant_name: str = "1a",
    product_name: str = "2a",
    write_results: bool = False,
    log_prefix: str = "",
    handle_errors: bool = False,
) -> None:
    if local_minima_path is None:
        local_minima_path = os.path.join(output_path, "local_minima")
    os.makedirs(os.path.join(local_minima_path, "reactant"), exist_ok=True)
    os.makedirs(os.path.join(local_minima_path, "product"), exist_ok=True)

    reactant_names = [reactant_name]
    product_names = [product_name]
    current_reactant_name = reactant_name
    current_product_name = product_name
    current_reactant_path = initial_reactant_path
    last_product_path = None

    csv_path = os.path.join(output_path, "scan.csv")
    csv_fp = open(csv_path, "w")
    csv_fp.write("reactant,product\n")

    try:
        while True:
            reaction_name = f"{current_reactant_name}-{current_product_name}"
            elementary_reaction_path = os.path.join(output_path, reaction_name)
            csv_fp.write(f"{current_reactant_name},{current_product_name}\n")
            csv_fp.flush()

            elementary_reaction_info = run_elementary_reaction(
                current_reactant_path,
                elementary_reaction_path,
                target_structure_path=last_product_path,
            )
            copy_reaction_local_minima(
                local_minima_path,
                current_reactant_name,
                current_product_name,
                elementary_reaction_path,
            )
            last_product_path = os.path.join(elementary_reaction_path, "product.xyz")

            intermediate_indices = elementary_reaction_info["intermediate_indices"]
            mep_path_info = elementary_reaction_info["mep_path_info"]
            ci_index = elementary_reaction_info["ci_index"]

            if elementary_reaction_info["terminate_chain"]:
                logger.info(
                    f"{log_prefix}CI at scan start (index 0) for {reaction_name} "
                    "(scan energy decreases along the path); no further scan needed"
                )
                if write_results:
                    write_rate_determining_ts_results(
                        output_path,
                        reaction_name,
                        mep_path_info,
                        ci_index,
                        local_minima_path,
                    )
                break

            chain_reactant_index = elementary_reaction_info["chain_reactant_index"]
            if chain_reactant_index is not None:
                logger.info(f"{log_prefix}Intermediate indices in the chain: {intermediate_indices}")
                current_reactant_path = os.path.join(
                    elementary_reaction_path,
                    f"{reaction_name}-{chain_reactant_index}.xyz",
                )
                ase.io.write(
                    current_reactant_path,
                    mep_path_info["atoms"][chain_reactant_index],
                    format="extxyz",
                )
                current_reactant_name = find_new_name(reactant_names)
                logger.info(
                    f"{log_prefix}New reactant {current_reactant_name} from interior minimum "
                    f"image {chain_reactant_index} of reaction {reaction_name} (CI at {ci_index})"
                )
                current_product_name = find_new_name(product_names)
                reactant_names.append(current_reactant_name)
                product_names.append(current_product_name)
            else:
                # No reactant-side interior minimum to chain from. Product-side
                # minima (if any) already define the product well via product_index;
                # treat the elementary scan as converged and keep the CI as RDTS.
                if intermediate_indices:
                    logger.info(
                        f"{log_prefix}Scan converged for reaction {reaction_name} "
                        f"(product-side minima at {intermediate_indices}; CI at {ci_index}; "
                        "no reactant-side minimum to chain)"
                    )
                else:
                    logger.info(f"{log_prefix}Scan converged for reaction {reaction_name}")
                if write_results:
                    write_rate_determining_ts_results(
                        output_path,
                        reaction_name,
                        mep_path_info,
                        ci_index,
                        local_minima_path,
                    )
                break
    except Exception:
        csv_fp.close()
        if handle_errors:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback_strs = traceback.format_exception(
                exc_type, exc_value, exc_traceback
            )
            logger.error(f"Error: {exc_type.__name__}: {exc_value}")
            logger.error(f"Error traceback: {''.join(traceback_strs)}")
            return
        raise
    finally:
        if not csv_fp.closed:
            csv_fp.close()
