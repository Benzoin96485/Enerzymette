import subprocess, sys, traceback, os, shutil
from typing import Optional
import yaml, json
from glob import glob
import ase
import ase.io
from ase.units import kcal, mol
from ..logger import logger
from ..mep_util import find_intermediate_indices, find_ci_index, find_rate_determining_step, find_new_name

class EnerzymeScanLauncher:
    def __init__(self,
        reactant_path: str,
        output_path: str,
        model_path: str,
        reference_path: str,
        reactant_name: str="1a",
        product_name: str="2a",
        n_steps: int=25,
    ):
        self.reactant_path = reactant_path
        self.output_path = output_path
        self.local_minima_path = os.path.join(self.output_path, "local_minima")
        self.model_path = model_path
        self.reference_path = reference_path
        self.reactant_name = reactant_name
        self.product_name = product_name
        self.n_steps = n_steps
        self.reference = self.parse_reference(reference_path, "terachem_input")
        self.constraint_freeze_xyz = self.reference.get("constraint_freeze", {}).get("xyz", [])
        self.constraint_scan = self.reference.get("constraint_scan", {})
        logger.info(f"Constraint freeze xyz: {self.constraint_freeze_xyz}")
        self.charge = int(self.reference.get("main", {}).get("charge", 0))
        logger.info(f"Charge: {self.charge}")
        self.multiplicity = int(self.reference.get("main", {}).get("spinmult", 1))
        logger.info(f"Multiplicity: {self.multiplicity}")

    def parse_reference(self, reference_path: str, reference_type: str):
        if reference_type == "terachem_input":
            logger.info(f"Parsing system information from terachem input file {reference_path}")
            from ..terachem.io import parse_terachem_input
            if not os.path.exists(reference_path):
                logger.warning(f"Reference file {reference_path} does not exist")
            else:
                return parse_terachem_input(reference_path)
        else:
            raise NotImplementedError(f"Reference type {reference_type} is not supported")
    
    def copy_local_minima(self, reactant_name: str, product_name: str):
        '''
        Copy the reactant and product of a scanned reaction with given reactant and product names to the local minima directory.
        '''
        elementary_reaction_path = os.path.join(self.output_path, f"{reactant_name}-{product_name}")
        reactant_path = os.path.join(elementary_reaction_path, f"reactant.xyz")
        product_path = os.path.join(elementary_reaction_path, f"product.xyz")
        if os.path.exists(reactant_path):
            shutil.copy(reactant_path, os.path.join(self.local_minima_path, "reactant", f"{reactant_name}.xyz"))
        else:
            logger.warning(f"Reactant file does not exist for finished reaction {reactant_name}-{product_name}")
        if os.path.exists(product_path):
            shutil.copy(product_path, os.path.join(self.local_minima_path, "product", f"{product_name}.xyz"))
        else:
            logger.warning(f"Product file does not exist for finished reaction {reactant_name}-{product_name}")

    def find_lowest_local_minima(self):
        '''
        Find the reactant and the product with the lowest energy from all reactant and product local minima of the entire reaction path.
        '''
        reactant_paths = glob(os.path.join(self.local_minima_path, "reactant", "*.xyz"))
        product_paths = glob(os.path.join(self.local_minima_path, "product", "*.xyz"))
        reactant_energies = [ase.io.read(reactant_path, format="extxyz", index=-1).get_potential_energy() for reactant_path in reactant_paths]
        product_energies = [ase.io.read(product_path, format="extxyz", index=-1).get_potential_energy() for product_path in product_paths]
        lowest_reactant_energy = min(reactant_energies)
        lowest_product_energy = min(product_energies)
        lowest_reactant_name = os.path.basename(reactant_paths[reactant_energies.index(lowest_reactant_energy)]).split(".")[0]
        lowest_product_name = os.path.basename(product_paths[product_energies.index(lowest_product_energy)]).split(".")[0]
        return lowest_reactant_name, lowest_product_name, lowest_reactant_energy, lowest_product_energy

    def launch(self):
        os.makedirs(self.local_minima_path, exist_ok=True)
        os.makedirs(os.path.join(self.local_minima_path, "reactant"), exist_ok=True)
        os.makedirs(os.path.join(self.local_minima_path, "product"), exist_ok=True)

        reactant_names = [self.reactant_name]
        product_names = [self.product_name]
        current_reactant_name = self.reactant_name
        current_product_name = self.product_name
        current_reactant_path = self.reactant_path
        last_product_path = None

        csv_path = os.path.join(self.output_path, "scan.csv")
        csv_fp = open(csv_path, "w")
        csv_fp.write("reactant,product\n")
        csv_fp.flush()

        try:
            while True:
                reaction_name = f"{current_reactant_name}-{current_product_name}"
                elementary_reaction_path = os.path.join(self.output_path, reaction_name)
                csv_fp.write(f"{current_reactant_name},{current_product_name}\n")
                csv_fp.flush()

                elementary_reaction_info = self.launch_elementary_reaction(current_reactant_path, elementary_reaction_path, target_structure_path=last_product_path)
                self.copy_local_minima(current_reactant_name, current_product_name)
                last_product_path = os.path.join(elementary_reaction_path, "product.xyz")
                
                intermediate_indices = elementary_reaction_info["intermediate_indices"]
                mep_path_info = elementary_reaction_info["mep_path_info"]
                ci_index = elementary_reaction_info["ci_index"]

                if len(intermediate_indices) > 0:
                    logger.info(f"Intermediate indices in the chain: {intermediate_indices}")
                    reactant_index, _, _ = find_rate_determining_step(intermediate_indices, mep_path_info["energies"], ci_index)
                    current_reactant_path = os.path.join(elementary_reaction_path, f"{reaction_name}-{reactant_index}.xyz")
                    ase.io.write(current_reactant_path, mep_path_info["atoms"][reactant_index], format="extxyz")
                    if reactant_index != 0:
                        current_reactant_name = find_new_name(reactant_names)
                        logger.info(f"New reactant {current_reactant_name} of rate determining step is found from image {reactant_index} of reaction {reaction_name}")

                    current_product_name = find_new_name(product_names)
                        
                    reactant_names.append(current_reactant_name)
                    product_names.append(current_product_name)
                else:
                    logger.info(f"Scan converged for reaction {reaction_name}")
                    lowest_reactant_name, lowest_product_name, lowest_reactant_energy, lowest_product_energy = self.find_lowest_local_minima()
                    ts_energy = mep_path_info["energies"][ci_index]
                    energy_span = (ts_energy - lowest_reactant_energy) / (kcal / mol)
                    energy_change = (lowest_product_energy - lowest_reactant_energy) / (kcal / mol)
                    logger.info(f"Reaction energy span: {energy_span:.2f} kcal/mol between reactant {lowest_reactant_name} and transition state of reaction {reaction_name}")
                    logger.info(f"Reaction energy change: {energy_change:.2f} kcal/mol between reactant {lowest_reactant_name} and product {lowest_product_name}")
                    results = {
                        "energy span": energy_span,
                        "energy change": energy_change,
                        "lowest energy reactant": lowest_reactant_name,
                        "lowest energy product": lowest_product_name
                    }
                    ts_dir = os.path.join(self.output_path, "rate_determining_ts")
                    os.makedirs(ts_dir, exist_ok=True)
                    with open(os.path.join(ts_dir, "results.json"), "w") as f:
                        json.dump(results, f, indent=4)
                    ase.io.write(os.path.join(ts_dir, f"{reaction_name}.xyz"), mep_path_info["atoms"][ci_index], format="extxyz")
                    break
        
        except Exception as e:
            csv_fp.close()
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback_strs = traceback.format_exception(exc_type, exc_value, exc_traceback)
            logger.error(f"Error: {exc_type.__name__}: {exc_value}")
            logger.error(f"Error traceback: {"".join(traceback_strs)}")

    def launch_elementary_reaction(self,
        reactant_path: str,
        elementary_reaction_path: str,
        target_value: Optional[float]=None,
        target_structure_path: Optional[str]=None,
    ):
        if os.path.exists(elementary_reaction_path):
            logger.warning(f"Output path {elementary_reaction_path} already exists, it could be overwritten")
        else:
            os.makedirs(elementary_reaction_path)

        logger.info(f"Reading reactant from {reactant_path}")
        reactant_atoms = ase.io.read(reactant_path, index=-1)
        init_reactant_path = os.path.join(elementary_reaction_path, "init_reactant.xyz")
        ase.io.write(init_reactant_path, reactant_atoms, format="extxyz")
        logger.info(f"Initial reactant written to {init_reactant_path}")

        # optimize reactant
        reactant_opt_config_path = os.path.join(elementary_reaction_path, "reactant_opt.yaml")
        self.write_config(
            task="opt",
            initial_structure_path=init_reactant_path,
            config_path=reactant_opt_config_path
        )
        opt_subprocess = subprocess.Popen(
            ["enerzyme", "simulate", "-c", reactant_opt_config_path, "-o", elementary_reaction_path, "-m", self.model_path]
        )
        opt_subprocess.wait()
        opt_path = os.path.join(elementary_reaction_path, "optim.xyz")
        reactant_path = os.path.join(elementary_reaction_path, "reactant.xyz")
        os.rename(opt_path, reactant_path)
        logger.info(f"Reactant optimized and written to {reactant_path}")
        
        # flexible scan
        scan_config_path = os.path.join(elementary_reaction_path, "scan.yaml")
        self.write_config(
            task="scan",
            initial_structure_path=reactant_path,
            config_path=scan_config_path,
            target_value=target_value,
            target_structure_path=target_structure_path
        )
        scan_subprocess = subprocess.Popen(
            ["enerzyme", "simulate", "-c", scan_config_path, "-o", elementary_reaction_path, "-m", self.model_path]
        )
        scan_subprocess.wait()
        logger.info(f"Scan finished for elementary reaction {elementary_reaction_path}")
        scan_atoms = ase.io.read(os.path.join(elementary_reaction_path, "scan_optim.xyz"), format="extxyz", index=":")
        mep_path_info = {
            "atoms": scan_atoms,
            "energies": [atoms.get_potential_energy() for atoms in scan_atoms],
            "n_images": len(scan_atoms),
            "n_atoms": len(scan_atoms[0])
        }

        intermediate_indices = find_intermediate_indices(mep_path_info["energies"])
        ci_index = find_ci_index(mep_path_info["energies"])
        _, product_index, _ = find_rate_determining_step(intermediate_indices, mep_path_info["energies"], ci_index)

        # optimize product
        scan_product_path = os.path.join(elementary_reaction_path, "init_product.xyz")
        ase.io.write(scan_product_path, scan_atoms[product_index], format="extxyz")
        logger.info(f"Scanned product (image {product_index}) written to {scan_product_path}")
        product_opt_config_path = os.path.join(elementary_reaction_path, "product_opt.yaml")
        self.write_config(
            task="opt",
            initial_structure_path=scan_product_path,
            config_path=product_opt_config_path
        )
        product_subprocess = subprocess.Popen(
            ["enerzyme", "simulate", "-c", product_opt_config_path, "-o", elementary_reaction_path, "-m", self.model_path],
        )
        product_subprocess.wait()
        product_path = os.path.join(elementary_reaction_path, "product.xyz")
        os.rename(opt_path, product_path)
        logger.info(f"Product optimized and written to {product_path}")

        return {
            "intermediate_indices": intermediate_indices,
            "mep_path_info": mep_path_info,
            "ci_index": ci_index,
        }

    def write_config(self, 
        task: str,
        initial_structure_path: str,
        config_path: str,
        target_value: Optional[float]=None,
        target_structure_path: Optional[str]=None,
    ):
        base_config = {
            "Simulation": {
                "environment": "ase",
                "dtype": "float64",
                "cuda": True,
                "task": task,
                "idx_start_from": 1,
                "neighbor_list": "full",
                "constraint": {
                    "fix_atom": {
                        "indices": self.constraint_freeze_xyz
                    }
                },
                "optimize": {
                    "optimizer": "LBFGS"
                }
            },
            "System": {
                "structure_file": initial_structure_path,
                "charge": self.charge,
                "multiplicity": self.multiplicity,
            }
        }
        if task == "scan":
            for constraint_type in self.constraint_scan.keys():
                constraint_params = self.constraint_scan[constraint_type]
                if constraint_type == "bond":
                    index0 = constraint_params["i0"]
                    index1 = constraint_params["i1"]
                    initial_structure = ase.io.read(initial_structure_path, index=-1)
                    initial_value = initial_structure.get_distance(index0 - 1, index1 - 1)
                    if target_value is None:
                        if target_structure_path is None:
                            from mendeleev import element
                            element0 = element(initial_structure.symbols[index0 - 1])
                            element1 = element(initial_structure.symbols[index1 - 1])
                            target_value = (element0.covalent_radius_pyykko + element1.covalent_radius_pyykko) / 100 # (pm to Angstrom)
                            logger.info(f"Target value is set to {target_value} Angstrom based on single-bond Pyykko covalent radii for bond between atoms {index0} and {index1}")
                        else:
                            target_structure = ase.io.read(target_structure_path, index=-1)
                            target_value = target_structure.get_distance(index0 - 1, index1 - 1)
                            logger.info(f"Target value is set to {target_value} Angstrom based on distance between atoms {index0} and {index1} in reference target structure {target_structure_path}")
                    base_config["Simulation"]["sampling"] = {
                        "cv": "distance",
                        "params": {
                            "x0": float(initial_value),
                            "x1": float(target_value),
                            "num": self.n_steps,
                            "i0": index0,
                            "i1": index1
                        }
                    }
                    logger.info(f"Scanning distance between atoms {index0} and {index1} from {initial_value} to {target_value} with {self.n_steps} steps")

        with open(config_path, "w") as f:
            yaml.dump(base_config, f)
        logger.info(f"Config (task: {task}) written to {config_path}")
