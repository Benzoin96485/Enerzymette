import os, subprocess, pickle, copy, sys, shutil, glob, json
from random import shuffle
from typing import Tuple, Optional, Dict, Callable, Literal, List
import yaml
import ase.io
from ..external_calculator import get_calculator_patch
from ..plumed_config_generator import (
    get_plumed_patch,
    get_scan_config_generator_name,
    resolve_scan_endpoints,
)
from ..mep_util import (
    find_intermediate_indices,
    find_ci_index,
    find_rate_determining_step,
    find_new_name,
)
from ..logger import logger


def update_config(old_config: Dict, new_config: Dict):
    for key, value in new_config.items():
        if key in old_config:
            if isinstance(value, dict) and isinstance(old_config[key], dict):
                update_config(old_config[key], value)
            else:
                old_config[key] = value
        else:
            old_config[key] = value


def collect_trajectory(simulation_trajectory_path: str, collection_path: str) -> None:
    traj = ase.io.read(simulation_trajectory_path, format="extxyz", index="1:")
    datapoints = []
    for frame in traj:
        datapoints.append(
            {
                "Za": frame.get_atomic_numbers(),
                "Ra": frame.get_positions(),
                "Q": frame.info.get("charge", 0),
                "S": frame.info.get("spin", 1) - 1
            }
        )
    with open(collection_path, "wb") as f:
        pickle.dump(datapoints, f)


resname_list = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
    "MSE": "M",
    "HID": "H"
}.keys()

class active_learning_launcher:
    def __init__(self,
        output_path: str,
        tmp_path: str,
        calculator_patch_key: str,
        plumed_patch_key: str,
        simulation_config_path: str,
        extraction_config_path: str,
        annotation_config_path: str,
        training_config_path: str,
        n_iterations: int,
        pretrain_path: Optional[str]=None,
        model_config_path: Optional[str]=None,
        reference_pdb_path: Optional[str]=None,
        template_sdf_path: Optional[str]=None,
        initial_xyz_path: Optional[str]=None,
        restraint_mode: Literal["hard", "soft"]="hard",
        training_ratio: float=0.8,
        cluster_inference_batch_size: int=4,
        n_presimulation_steps_per_iteration: int=0,
        continual_learning: bool=False,
        reset_parameters: bool=False,
        initial_scan: bool=False,
        n_initial_scan_steps: int=25,
    ) -> None:
        self.pretrain_path = pretrain_path
        self.pretrain_config = None
        self.model_config_path = model_config_path
        self.reset_parameters = reset_parameters
        self.output_path = os.path.abspath(output_path)
        self.tmp_path = tmp_path
        self.calculator_patch = get_calculator_patch(calculator_patch_key)
        self.plumed_patch_key = plumed_patch_key
        self.plumed_patch = get_plumed_patch(plumed_patch_key)
        self.simulation_config_path = simulation_config_path
        self.extraction_config_path = extraction_config_path
        self.annotation_config_path = annotation_config_path
        self.training_config_path = training_config_path
        self.n_iterations = n_iterations
        self.model_architecture = None
        self.model_suffix = None
        self.active_model_key = None
        self.model_str = None
        self.training_ratio = training_ratio
        self.cluster_inference_batch_size = cluster_inference_batch_size
        self.n_presimulation_steps_per_iteration = n_presimulation_steps_per_iteration
        self.initial_scan = initial_scan
        self.n_initial_scan_steps = n_initial_scan_steps
        self.structure_pool_state: Optional[Dict] = None
        self.NA = None
        self.continual_learning = continual_learning
        self.reference_pdb_path = reference_pdb_path
        self.template_sdf_path = template_sdf_path
        self.restraint_mode = restraint_mode
        self.charge = None
        self.backbone_indices = None
        self.Calpha_indices = None
        self.output_img_path = os.path.join(self.output_path, "cluster.png")
        self.output_xyz_path = os.path.join(self.output_path, "cluster.xyz")
        self.output_mol_path = os.path.join(self.output_path, "cluster.mol")
        self.initial_xyz_path = initial_xyz_path

    def _get_topology(self) -> None:
        if self.reference_pdb_path is not None:
            from rdkit.Chem import MolFromMolFile, MolToXYZFile, GetFormalCharge
            enerzyme_bond_args = [
                "enerzyme", "bond",
                "-p", self.reference_pdb_path,
                "-i", self.output_img_path,
                "-m", self.output_mol_path,
            ]
            if self.template_sdf_path is not None:
                enerzyme_bond_args.extend(["-t", self.template_sdf_path])
            enerzyme_subprocess = subprocess.Popen(enerzyme_bond_args)
            enerzyme_subprocess.wait()
            cluster_mol = MolFromMolFile(self.output_mol_path, removeHs=False)
            self.charge = GetFormalCharge(cluster_mol)
            if self.initial_xyz_path is None:
                MolToXYZFile(cluster_mol, self.output_xyz_path)
                self.initial_xyz_path = self.output_xyz_path

            self.backbone_indices = []
            self.Calpha_indices = []
            atom_count = 0
            with open(self.reference_pdb_path, "r") as f:
                for line in f.readlines():
                    if line.startswith("ATOM") or line.startswith("HETATM"):
                        atom_count += 1
                        resname = line[17:20]
                        atomname = line[11:16].strip()
                        if resname in resname_list:
                            if atomname in ["N", "C", "CA", "O"]:
                                self.backbone_indices.append(atom_count - 1)
                                if atomname == "CA":
                                    self.Calpha_indices.append(atom_count - 1)
                        elif resname == "HOH":
                            if atomname == "O":
                                self.Calpha_indices.append(atom_count - 1)
            with open(self.simulation_config_path, "r") as f:
                simulation_config = yaml.load(f, Loader=yaml.FullLoader)
                simulation_config["System"]["structure_file"] = self.initial_xyz_path
                simulation_config["System"]["charge"] = self.charge
                idx_start_from = simulation_config["Simulation"]["idx_start_from"]
                if self.restraint_mode == "hard":
                    simulation_config["Simulation"]["constraint"] = {
                        "fix_atom": {
                            "indices": [idx + idx_start_from for idx in self.backbone_indices]
                        }
                    }
                elif self.restraint_mode == "soft":
                    Hookean_k = simulation_config["Simulation"]["constraint"].get("Hookean_allpairs", {}).get("k", 0.05)
                    simulation_config["Simulation"]["constraint"] = {
                        "Hookean_allpairs": {
                            "indices": [idx + idx_start_from for idx in self.Calpha_indices],
                            "k": Hookean_k
                        }
                    }
                simulation_config["Simulation"]["sampling"]["params"]["plumed_config"]["reference_pdb_file"] = self.reference_pdb_path
            with open(self.simulation_config_path, "w") as f:
                yaml.dump(simulation_config, f, default_flow_style=False)
            logger.info(f"Using total charge ({self.charge}) of the cluster parsed from the reference structure: {self.reference_pdb_path} with the small molecule template: {self.template_sdf_path}")
            logger.info(f"Using initial structure: {self.initial_xyz_path} converted from the reference structure: {self.reference_pdb_path}")
            
            with open(self.extraction_config_path, "r") as f:
                extraction_config = yaml.load(f, Loader=yaml.FullLoader)
                extraction_config["Extractor"]["reference_mol_path"] = self.output_mol_path
            with open(self.extraction_config_path, "w") as f:
                yaml.dump(extraction_config, f, default_flow_style=False)
            logger.info(f"Using topology information {self.output_mol_path} parsed from the reference structure: {self.reference_pdb_path} and the small molecule template: {self.template_sdf_path}")


    def _get_model_suffix(self, i: int) -> str:
        if not self.model_suffix:
            return f"{i}"
        else:
            return f"{self.model_suffix}-{i}"

    def _get_model_config(self) -> None:
        if self.pretrain_path is not None:  
            pretrain_config_path = os.path.join(self.pretrain_path, "config.yaml")
        else:
            pretrain_config_path = self.model_config_path
        with open(pretrain_config_path, "r") as f:
            self.pretrain_config = yaml.load(f, Loader=yaml.FullLoader)
        active_model_key = None
        all_internal_FFs = self.pretrain_config["Modelhub"]["internal_FFs"]
        for model_key, model_params in all_internal_FFs.items():
            if model_params.get("active", False) == True:
                active_model_key = model_key
        if active_model_key is None:
            raise ValueError(f"No active model found in {pretrain_config_path}")
        else:
            self.active_model_key = active_model_key

        self.model_architecture = all_internal_FFs[self.active_model_key]["architecture"]
        self.model_suffix = all_internal_FFs[self.active_model_key].get("suffix", "")

    def _make_model_config(self, i: int, dump_path: str, 
        new_config_path: Optional[str]=None, 
        new_datahub_config: Optional[Dict]=None, 
        new_trainer_config: Optional[Dict]=None, 
        model_params_updater: Optional[Callable[[int, Dict], None]]=None,
        new_metric_config: Optional[Dict]=None,
    ) -> None:
        suffix = self._get_model_suffix(i)
        model_config = copy.deepcopy(self.pretrain_config)
        model_param = model_config["Modelhub"]["internal_FFs"][self.active_model_key]    
        if new_config_path is not None:
            with open(new_config_path, "r") as f:
                new_config = yaml.load(f, Loader=yaml.FullLoader)
            update_config(model_config, new_config)
        
        if new_datahub_config is not None:
            if model_config.get("Datahub", None) is None:
                model_config["Datahub"] = {}
            update_config(model_config["Datahub"], new_datahub_config)
        
        if new_trainer_config is not None:
            if model_config.get("Trainer", None) is None:
                model_config["Trainer"] = {}
            update_config(model_config["Trainer"], new_trainer_config)

        if model_params_updater is not None:
            model_params_updater(i, model_param)

        if new_metric_config is not None:
            if model_config.get("Metric", None) is None:
                model_config["Metric"] = {}
            model_config["Metric"] = copy.deepcopy(new_metric_config)

        model_param["suffix"] = suffix

        with open(dump_path, "w") as f:
            yaml.dump(model_config, f, default_flow_style=False)

    def _get_model_str(self, i: Optional[int]=None) -> str:
        if i is None:
            return f"{self.active_model_key}-{self.model_architecture}"
        else:
            return f"{self.active_model_key}-{self.model_architecture}-{self._get_model_suffix(i)}"

    def _get_prediction_path(self, i: int) -> Tuple[str, str, str, str]:
        prediction_path = os.path.join(self.output_path, self._get_model_str(i) + "_prediction")
        prediction_config_path = os.path.join(prediction_path, "prediction.yaml")
        prediction_model_config_path = os.path.join(prediction_path, "model_config.yaml")
        prediction_glob_pattern = os.path.join(prediction_path, f"processed_dataset_*/{self._get_model_str(i)}-prediction.pkl")
        return prediction_path, prediction_config_path, prediction_model_config_path, prediction_glob_pattern

    def _prediction_step(self, i: int) -> None:
        if i <= 0:
            logger.info(f"Direct simulation for iteration {i}")
            return
        prediction_path, prediction_config_path, prediction_model_config_path, prediction_glob_pattern = self._get_prediction_path(i)

        prediction_pickles = glob.glob(prediction_glob_pattern)
        if prediction_pickles:
            logger.info(f"Prediction found at {prediction_pickles[0]} and skipped.")
            return
        
        os.makedirs(prediction_path, exist_ok=True)

        collection_path = self._get_collection_path(i - 1)
        prediction_config = {
            "Datahub": copy.deepcopy(self.pretrain_config["Datahub"]),
        }
        prediction_config["Datahub"].pop("targets")
        prediction_config["Datahub"]["data_path"] = collection_path
        prediction_config["Datahub"]["features"] = {
            "N": None,
            "Za": "Za",
            "Ra": "Ra",
            "Q": "Q"
        }
        with open(prediction_config_path, "w") as f:
            yaml.dump(prediction_config, f, default_flow_style=False)

        self._make_model_config(i, prediction_model_config_path, new_trainer_config={"inference_batch_size": self.cluster_inference_batch_size})

        enerzyme_subprocess = subprocess.Popen(
            ["enerzyme", "predict", 
                "-c", prediction_config_path, 
                "-o", prediction_path, 
                "-m", self.output_path,
                "-mc", prediction_model_config_path,
                "-s"
            ],
            stdout=sys.stdout,
            stderr=sys.stderr
        )
        enerzyme_subprocess.wait()
        logger.info(f"Prediction finished for {prediction_path}")

    def _get_simulation_path(self, i: int) -> Tuple[str, str, str, str, str]:
        simulation_path = os.path.join(self.output_path, f"{self._get_model_str(i)}_md")
        simulation_trajectory_path = os.path.join(simulation_path, "plumed.traj.xyz")
        simulation_config_path = os.path.join(simulation_path, "simulation.yaml")
        simulation_model_config_path = os.path.join(simulation_path, "model_config.yaml")
        simulation_completed_flag = os.path.join(simulation_path, "simulation_completed")
        initial_structure_path = os.path.join(simulation_path, "initial_structure.xyz")
        presimulation_config_path = os.path.join(simulation_path,
            "presimulation.yaml"
        )
        presimulation_trajectory_path = os.path.join(simulation_path, "md.traj.xyz")
        presimulation_completed_flag = os.path.join(simulation_path, "presimulation_completed")
        return simulation_path, simulation_trajectory_path, simulation_config_path, simulation_model_config_path, simulation_completed_flag, initial_structure_path, presimulation_config_path, presimulation_trajectory_path, presimulation_completed_flag

    def _get_max_E_var(self, i: int) -> float:
        _, _, _, prediction_glob_pattern = self._get_prediction_path(i)
        prediction_pickles = glob.glob(prediction_glob_pattern)
        if prediction_pickles:
            prediction_pickle = prediction_pickles[0]
            with open(prediction_pickle, "rb") as f:
                prediction_data = pickle.load(f)
            print(prediction_data.keys())
            average_E_var = prediction_data["E_var"].mean().item()
            logger.info(f"Average E_var for {prediction_pickle}: {average_E_var}")
            return average_E_var
        else:
            raise FileNotFoundError(f"Prediction pickles not found for {prediction_glob_pattern}")
    
    def _udd_model_params_updater(self, i: int, model_param: Dict) -> None:
        shallow_ensemble_reduce_found = False
        new_layers = []
        for layer in model_param["layers"]:
            if layer["name"] == "ShallowEnsembleReduce":
                if not shallow_ensemble_reduce_found:
                    shallow_ensemble_reduce_found = True
                    if "train_only" in layer["params"]:
                        layer["params"].pop("train_only")
                else:
                    continue
            elif layer["name"] == "Force":
                new_layers.append({"name": "EnergyVarianceGradient"})
            new_layers.append(layer)
        model_param["layers"] = new_layers

    def _structure_pool_state_path(self) -> str:
        return os.path.join(self.output_path, "structure_pool.json")

    def _structure_pool_dir(self) -> str:
        return os.path.join(self.output_path, "structure_pool")

    def _load_structure_pool(self) -> Dict:
        state_path = self._structure_pool_state_path()
        if os.path.exists(state_path):
            with open(state_path, "r") as f:
                return json.load(f)
        return {"entries": []}

    def _save_structure_pool(self) -> None:
        with open(self._structure_pool_state_path(), "w") as f:
            json.dump(self.structure_pool_state, f, indent=2)

    def _write_task_config(
        self,
        task: str,
        structure_path: str,
        config_path: str,
        base_config: Dict,
        scan_target_value: Optional[float] = None,
        scan_target_structure_path: Optional[str] = None,
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
        elif task == "md":
            config["Simulation"].pop("sampling", None)
            config["Simulation"]["integrate"]["n_step"] = self.n_presimulation_steps_per_iteration
        elif task == "plumed_scan":
            config["Simulation"].pop("integrate", None)
            if "optimize" not in config["Simulation"]:
                config["Simulation"]["optimize"] = {"optimizer": "LBFGS"}
            idx_start_from = config["Simulation"]["idx_start_from"]
            plumed_cv_config = config["Simulation"]["sampling"]["params"]["plumed_config"]
            initial_structure = ase.io.read(structure_path, index=-1)
            x0, x1, num, rc = resolve_scan_endpoints(
                initial_structure,
                idx_start_from,
                self.plumed_patch_key,
                plumed_cv_config,
                self.n_initial_scan_steps,
                target_value=scan_target_value,
                target_structure_path=scan_target_structure_path,
            )
            config["Simulation"]["plumed_config_generator"] = {
                "name": get_scan_config_generator_name(self.plumed_patch_key),
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

        with open(config_path, "w") as f:
            yaml.dump(config, f, default_flow_style=False)

    def _enerzyme_simulate_cmd(
        self,
        config_path: str,
        output_path: str,
        model_config_path: str,
        with_plumed_patch: bool = False,
    ) -> List[str]:
        cmd = [
            "enerzyme", "simulate",
            "-c", config_path,
            "-o", output_path,
            "-m", self.output_path,
            "-cp", self.calculator_patch,
            "-mc", model_config_path,
        ]
        if with_plumed_patch:
            cmd.extend(["-pp", self.plumed_patch])
        return cmd

    def _copy_initial_scan_local_minima(
        self,
        local_minima_path: str,
        reactant_name: str,
        product_name: str,
        elementary_reaction_path: str,
    ) -> None:
        reactant_src = os.path.join(elementary_reaction_path, "reactant.xyz")
        product_src = os.path.join(elementary_reaction_path, "product.xyz")
        if os.path.exists(reactant_src):
            shutil.copy(
                reactant_src,
                os.path.join(local_minima_path, "reactant", f"{reactant_name}.xyz"),
            )
        if os.path.exists(product_src):
            shutil.copy(
                product_src,
                os.path.join(local_minima_path, "product", f"{product_name}.xyz"),
            )

    def _launch_elementary_reaction_scan(
        self,
        reactant_path: str,
        elementary_reaction_path: str,
        base_config: Dict,
        model_config_path: str,
        target_value: Optional[float] = None,
        target_structure_path: Optional[str] = None,
    ) -> Dict:
        if os.path.exists(elementary_reaction_path):
            logger.warning(
                f"Output path {elementary_reaction_path} already exists, it could be overwritten"
            )
        else:
            os.makedirs(elementary_reaction_path)

        reactant_atoms = ase.io.read(reactant_path, index=-1)
        init_reactant_path = os.path.join(elementary_reaction_path, "init_reactant.xyz")
        ase.io.write(init_reactant_path, reactant_atoms, format="extxyz")

        reactant_opt_config_path = os.path.join(elementary_reaction_path, "reactant_opt.yaml")
        self._write_task_config("opt", init_reactant_path, reactant_opt_config_path, base_config)
        subprocess.Popen(
            self._enerzyme_simulate_cmd(reactant_opt_config_path, elementary_reaction_path, model_config_path),
            stdout=sys.stdout,
            stderr=sys.stderr,
        ).wait()
        opt_path = os.path.join(elementary_reaction_path, "optim.xyz")
        optimized_reactant_path = os.path.join(elementary_reaction_path, "reactant.xyz")
        os.rename(opt_path, optimized_reactant_path)

        scan_config_path = os.path.join(elementary_reaction_path, "scan.yaml")
        self._write_task_config(
            "plumed_scan",
            optimized_reactant_path,
            scan_config_path,
            base_config,
            scan_target_value=target_value,
            scan_target_structure_path=target_structure_path,
        )
        subprocess.Popen(
            self._enerzyme_simulate_cmd(
                scan_config_path, elementary_reaction_path, model_config_path, with_plumed_patch=True
            ),
            stdout=sys.stdout,
            stderr=sys.stderr,
        ).wait()

        scan_atoms = ase.io.read(
            os.path.join(elementary_reaction_path, "scan_optim.xyz"), format="extxyz", index=":"
        )
        energies = [atoms.get_potential_energy() for atoms in scan_atoms]
        intermediate_indices = find_intermediate_indices(energies)
        ci_index = find_ci_index(energies)
        _, product_index, _ = find_rate_determining_step(intermediate_indices, energies, ci_index)

        scan_product_path = os.path.join(elementary_reaction_path, "init_product.xyz")
        ase.io.write(scan_product_path, scan_atoms[product_index], format="extxyz")
        product_opt_config_path = os.path.join(elementary_reaction_path, "product_opt.yaml")
        self._write_task_config("opt", scan_product_path, product_opt_config_path, base_config)
        subprocess.Popen(
            self._enerzyme_simulate_cmd(product_opt_config_path, elementary_reaction_path, model_config_path),
            stdout=sys.stdout,
            stderr=sys.stderr,
        ).wait()
        os.rename(opt_path, os.path.join(elementary_reaction_path, "product.xyz"))

        return {
            "intermediate_indices": intermediate_indices,
            "mep_path_info": {"atoms": scan_atoms, "energies": energies},
            "ci_index": ci_index,
        }

    def _run_initial_scan(self) -> None:
        initial_scan_path = os.path.join(self.output_path, "initial_scan")
        local_minima_path = os.path.join(initial_scan_path, "local_minima")
        completed_flag = os.path.join(initial_scan_path, "initial_scan_completed")
        if os.path.exists(completed_flag):
            logger.info(f"Initial scan already completed at {initial_scan_path}, skipping.")
            return

        os.makedirs(os.path.join(local_minima_path, "reactant"), exist_ok=True)
        os.makedirs(os.path.join(local_minima_path, "product"), exist_ok=True)

        model_config_path = os.path.join(initial_scan_path, "model_config.yaml")
        if not os.path.exists(model_config_path):
            self._make_model_config(0, model_config_path)

        with open(self.simulation_config_path, "r") as f:
            base_config = yaml.load(f, Loader=yaml.FullLoader)

        reactant_names = ["1a"]
        product_names = ["2a"]
        current_reactant_name = reactant_names[0]
        current_product_name = product_names[0]
        current_reactant_path = base_config["System"]["structure_file"]
        last_product_path = None

        csv_path = os.path.join(initial_scan_path, "scan.csv")
        with open(csv_path, "w") as csv_fp:
            csv_fp.write("reactant,product\n")
            while True:
                reaction_name = f"{current_reactant_name}-{current_product_name}"
                elementary_reaction_path = os.path.join(initial_scan_path, reaction_name)
                csv_fp.write(f"{current_reactant_name},{current_product_name}\n")
                csv_fp.flush()

                elementary_reaction_info = self._launch_elementary_reaction_scan(
                    current_reactant_path,
                    elementary_reaction_path,
                    base_config,
                    model_config_path,
                    target_structure_path=last_product_path,
                )
                self._copy_initial_scan_local_minima(
                    local_minima_path,
                    current_reactant_name,
                    current_product_name,
                    elementary_reaction_path,
                )
                last_product_path = os.path.join(elementary_reaction_path, "product.xyz")

                intermediate_indices = elementary_reaction_info["intermediate_indices"]
                mep_path_info = elementary_reaction_info["mep_path_info"]
                ci_index = elementary_reaction_info["ci_index"]

                if len(intermediate_indices) > 0:
                    reactant_index, _, _ = find_rate_determining_step(
                        intermediate_indices, mep_path_info["energies"], ci_index
                    )
                    current_reactant_path = os.path.join(
                        elementary_reaction_path, f"{reaction_name}-{reactant_index}.xyz"
                    )
                    ase.io.write(
                        current_reactant_path,
                        mep_path_info["atoms"][reactant_index],
                        format="extxyz",
                    )
                    if reactant_index != 0:
                        current_reactant_name = find_new_name(reactant_names)
                        reactant_names.append(current_reactant_name)
                    current_product_name = find_new_name(product_names)
                    product_names.append(current_product_name)
                else:
                    logger.info(f"Initial scan converged for reaction {reaction_name}")
                    break

        open(completed_flag, "w").close()
        logger.info(f"Initial scan finished at {initial_scan_path}")

    def _build_structure_pool_from_initial_scan(self) -> None:
        initial_scan_path = os.path.join(self.output_path, "initial_scan")
        local_minima_path = os.path.join(initial_scan_path, "local_minima")
        pool_dir = self._structure_pool_dir()
        os.makedirs(pool_dir, exist_ok=True)

        pool_entries = []
        csv_path = os.path.join(initial_scan_path, "scan.csv")
        with open(csv_path, "r") as f:
            lines = f.readlines()[1:]
        pool_idx = 0
        for line in lines:
            reactant_name, product_name = line.strip().split(",")
            for name, subdir in ((reactant_name, "reactant"), (product_name, "product")):
                src = os.path.join(local_minima_path, subdir, f"{name}.xyz")
                if not os.path.exists(src):
                    logger.warning(f"Local minimum not found for initial scan: {src}")
                    continue
                dst = os.path.join(pool_dir, f"{pool_idx:03d}.xyz")
                shutil.copy(src, dst)
                pool_entries.append({"path": dst, "run": False})
                pool_idx += 1

        if not pool_entries:
            raise RuntimeError(
                f"Initial scan produced no local minima under {local_minima_path}"
            )
        self.structure_pool_state = {"entries": pool_entries}
        self._save_structure_pool()
        logger.info(
            f"Structure pool initialized with {len(pool_entries)} entries from initial scan"
        )

    def _init_structure_pool(self) -> None:
        state_path = self._structure_pool_state_path()
        if os.path.exists(state_path):
            self.structure_pool_state = self._load_structure_pool()
            logger.info(
                f"Loaded structure pool with {len(self.structure_pool_state['entries'])} entries"
            )
            return

        if self.initial_scan:
            self._run_initial_scan()
            self._build_structure_pool_from_initial_scan()
            return

        with open(self.simulation_config_path, "r") as f:
            simulation_config = yaml.load(f, Loader=yaml.FullLoader)
        initial_path = os.path.abspath(simulation_config["System"]["structure_file"])
        pool_dir = self._structure_pool_dir()
        os.makedirs(pool_dir, exist_ok=True)
        pool_path = os.path.join(pool_dir, "000.xyz")
        if os.path.abspath(initial_path) != os.path.abspath(pool_path):
            shutil.copy(initial_path, pool_path)
        self.structure_pool_state = {"entries": [{"path": pool_path, "run": False}]}
        self._save_structure_pool()
        logger.info(f"Structure pool initialized with initial structure: {pool_path}")

    def _prepare_iteration_initial_structure(
        self,
        i: int,
        initial_structure_path: str,
        simulation_path: str,
        presimulation_config_path: str,
        presimulation_trajectory_path: str,
        presimulation_completed_flag: str,
        simulation_model_config_path: str,
        simulation_config: Dict,
    ) -> None:
        entries = self.structure_pool_state["entries"]
        pool_idx = i % len(entries)
        entry = entries[pool_idx]

        if not entry["run"]:
            shutil.copy(entry["path"], initial_structure_path)
            entry["run"] = True
            self._save_structure_pool()
            logger.info(
                f"Iteration {i}: using structure pool entry {pool_idx} directly "
                f"({entry['path']})"
            )
            return

        if self.n_presimulation_steps_per_iteration > 0:
            if os.path.exists(presimulation_completed_flag):
                logger.info(f"Presimulation completed for {simulation_path} and skipped.")
            else:
                presimulation_config = copy.deepcopy(simulation_config)
                presimulation_config["Simulation"].pop("sampling", None)
                presimulation_config["Simulation"]["task"] = "md"
                presimulation_config["Simulation"]["integrate"]["n_step"] = (
                    self.n_presimulation_steps_per_iteration
                )
                presimulation_config["System"]["structure_file"] = entry["path"]
                with open(presimulation_config_path, "w") as f:
                    yaml.dump(presimulation_config, f, default_flow_style=False)
                subprocess.Popen(
                    self._enerzyme_simulate_cmd(
                        presimulation_config_path,
                        simulation_path,
                        simulation_model_config_path,
                    ),
                    stdout=sys.stdout,
                    stderr=sys.stderr,
                ).wait()
                if os.path.exists(presimulation_trajectory_path):
                    open(presimulation_completed_flag, "w").close()
            last_frame = ase.io.read(presimulation_trajectory_path, index=-1)
            ase.io.write(entry["path"], last_frame, format="extxyz")
            ase.io.write(initial_structure_path, last_frame, format="extxyz")
            self._save_structure_pool()
            logger.info(
                f"Iteration {i}: presimulation from pool entry {pool_idx}, "
                f"updated entry at {entry['path']}"
            )
            return

        last_simulation_trajectory_path = self._get_simulation_path(i - 1)[1]
        last_frame = ase.io.read(last_simulation_trajectory_path, index=-1)
        ase.io.write(initial_structure_path, last_frame, format="extxyz")
        logger.info(
            f"Iteration {i}: reusing last steered-MD frame from iteration {i - 1}"
        )

    def _simulation_step(self, i: int) -> None:
        simulation_path, simulation_trajectory_path, simulation_config_path, simulation_model_config_path, simulation_completed_flag, initial_structure_path, presimulation_config_path, presimulation_trajectory_path, presimulation_completed_flag = self._get_simulation_path(i)

        if os.path.exists(simulation_completed_flag):
            logger.info(f"Simulation completed for {simulation_path} and skipped.")
            return
        else:
            if os.path.exists(simulation_path):
                shutil.rmtree(simulation_path)
        
        os.makedirs(simulation_path, exist_ok=True)

        with open(self.simulation_config_path, "r") as f:
            simulation_config = yaml.load(f, Loader=yaml.FullLoader)
        if i == 0:
            simulation_config["Simulation"].pop("uncertainty_calculator", None)
        else:
            average_E_var = self._get_max_E_var(i)
            if self.NA is None:
                initial_structure = ase.io.read(simulation_config["System"]["structure_file"])
                self.NA = len(initial_structure)
                logger.info(f"Initial structure has {self.NA} atoms.")
            simulation_config["Simulation"]["uncertainty_calculator"]["params"]["B"] = average_E_var / self.NA

        simulation_config["System"]["structure_file"] = initial_structure_path

        with open(simulation_config_path, "w") as f:
            yaml.dump(simulation_config, f, default_flow_style=False)

        self._make_model_config(i, simulation_model_config_path, model_params_updater=self._udd_model_params_updater)

        self._prepare_iteration_initial_structure(
            i,
            initial_structure_path,
            simulation_path,
            presimulation_config_path,
            presimulation_trajectory_path,
            presimulation_completed_flag,
            simulation_model_config_path,
            simulation_config,
        )

        enerzyme_subprocess = subprocess.Popen(
            ["enerzyme", "simulate", 
                "-c", simulation_config_path, 
                "-o", simulation_path, 
                "-m", self.output_path, 
                "-cp", self.calculator_patch,
                "-pp", self.plumed_patch,
                "-mc", simulation_model_config_path
            ],
            stdout=sys.stdout,
            stderr=sys.stderr
        )
        enerzyme_subprocess.wait()

        if os.path.exists(simulation_trajectory_path):
            open(simulation_completed_flag, "w").close()
        logger.info(f"Simulation finished for {simulation_path}")

    def _get_collection_path(self, i: int) -> str:
        simulation_path, *_ = self._get_simulation_path(i)
        return os.path.join(simulation_path, "plumed.traj.pkl")

    def _collection_step(self, i: int) -> None:
        _, simulation_trajectory_path, *_ = self._get_simulation_path(i)
        collection_path = self._get_collection_path(i)
        if os.path.exists(collection_path):
            logger.info(f"Raw dataset found at {collection_path}. Collection skipped")
            return
        if os.path.exists(simulation_trajectory_path):
            collect_trajectory(simulation_trajectory_path, collection_path)
            logger.info(f"Raw dataset collected to {collection_path}")
        else:
            raise FileNotFoundError(f"Simulation trajectory not found for {simulation_trajectory_path}")

    def _get_extraction_path(self, i: int) -> Tuple[str, str, str, str]:
        extraction_path = os.path.join(self.output_path, self._get_model_str(i) + "_extraction")
        extraction_file_path = os.path.join(extraction_path, self._get_model_str(i) + "_fragments.sdf")
        extraction_config_path = os.path.join(extraction_path, "extraction.yaml")
        extraction_model_config_path = os.path.join(extraction_path, "model_config.yaml")
        return extraction_path, extraction_file_path, extraction_config_path, extraction_model_config_path

    def _extraction_step(self, i: int) -> None:
        extraction_path, extraction_file_path, extraction_config_path, extraction_model_config_path = self._get_extraction_path(i)
        if os.path.exists(extraction_file_path):
            logger.info(f"Fragments found at {extraction_file_path}. Extraction skipped")
            return
        
        os.makedirs(extraction_path, exist_ok=True)

        collection_path = self._get_collection_path(i)
        with open(self.extraction_config_path, "r") as f:
            extraction_config = yaml.load(f, Loader=yaml.FullLoader)
        extraction_config["Datahub"]["data_path"] = collection_path
        if i <= 0 and (self.reset_parameters or self.pretrain_path is None):
            extraction_config["Extractor"]["extract_method"] = "random"
        with open(extraction_config_path, "w") as f:
            yaml.dump(extraction_config, f, default_flow_style=False)

        new_trainer_config = {
            "inference_batch_size": self.cluster_inference_batch_size
        }
        self._make_model_config(i, extraction_model_config_path, new_trainer_config=new_trainer_config)

        enerzyme_subprocess = subprocess.Popen(
            ["enerzyme", "extract", 
                "-c", extraction_config_path, 
                "-o", extraction_path, 
                "-m", self.output_path,
                "-mc", extraction_model_config_path
            ],
            stdout=sys.stdout,
            stderr=sys.stderr
        )
        enerzyme_subprocess.wait()
        logger.info(f"Fragments extracted to {extraction_file_path}")

    def _get_annotation_path(self, i: int) -> str:
        extraction_path, *_ = self._get_extraction_path(i)
        annotation_path = os.path.join(self.output_path, self._get_model_str(i) + "_fragments")
        annotation_file_path = os.path.join(annotation_path, "fragments.pkl")
        annotation_config_path = os.path.join(extraction_path, "annotation.yaml")
        return annotation_file_path, annotation_config_path

    def _annotation_step(self, i: int) -> None:
        annotation_file_path, annotation_config_path = self._get_annotation_path(i)
        if os.path.exists(annotation_file_path):
            logger.info(f"Annotated dataset found at {annotation_file_path}. Annotation skipped.")
            return
        _, extraction_file_path, *_ = self._get_extraction_path(i)

        with open(self.annotation_config_path, "r") as f:
            annotation_config = yaml.load(f, Loader=yaml.FullLoader)
        annotation_config["Supplier"]["path"] = extraction_file_path
        with open(annotation_config_path, "w") as f:
            yaml.dump(annotation_config, f)
        
        enerzyme_subprocess = subprocess.Popen(
            ["enerzyme", "annotate", "-c", annotation_config_path, "-o", self.output_path, "-t", self.tmp_path],
            stdout=sys.stdout,
            stderr=sys.stderr
        )
        enerzyme_subprocess.wait()
        logger.info(f"Annotated dataset saved to {annotation_file_path}")

    def _get_training_path(self, i: int) -> str:
        training_path = os.path.join(self.output_path, self._get_model_str(i + 1) + "_training")
        training_config_path = os.path.join(training_path, "train.yaml")
        training_set_path = os.path.join(training_path, "training_set.pkl")
        validation_set_path = os.path.join(training_path, "validation_set.pkl")
        training_completed_flag = os.path.join(training_path, "training_completed")
        return training_path, training_config_path, training_set_path, validation_set_path, training_completed_flag

    def _get_merged_dataset_path(self, i: int) -> str:
        return os.path.join(self.output_path, f"merged_dataset_{i}.pkl")

    def _merge_dataset(self, i: int) -> None:
        _, _, training_set_path, validation_set_path, _ = self._get_training_path(i)
        if os.path.exists(training_set_path) and os.path.exists(validation_set_path):
            logger.info(f"Merged dataset found at {training_set_path} and {validation_set_path}. Merging skipped.")
            return
        
        if i > 0:
            _, _, old_training_set_path, old_validation_set_path, _ = self._get_training_path(i - 1)
            with open(old_training_set_path, "rb") as f:
                all_training_datapoints = pickle.load(f)
            with open(old_validation_set_path, "rb") as f:
                all_validation_datapoints = pickle.load(f)
        else:
            all_training_datapoints = []
            all_validation_datapoints = []
            
        annotation_file_path, _ = self._get_annotation_path(i)
        with open(annotation_file_path, "rb") as f:
            datapoints = pickle.load(f)
        shuffle(datapoints)
        len_training_set = int(len(datapoints) * self.training_ratio)
        all_training_datapoints.extend(datapoints[:len_training_set])
        all_validation_datapoints.extend(datapoints[len_training_set:])

        with open(training_set_path, "wb") as f:
            pickle.dump(all_training_datapoints, f)
        with open(validation_set_path, "wb") as f:
            pickle.dump(all_validation_datapoints, f)
        logger.info(f"Merged dataset saved to {training_set_path} and {validation_set_path}")

    def _training_model_params_updater(self, i: int, model_param: Dict) -> None:
        model_param["pretrain_path"] = os.path.join(self.output_path, self._get_model_str(i - 1))

    def _training_step(self, i: int) -> None:
        training_path, training_config_path, training_set_path, validation_set_path, training_completed_flag = self._get_training_path(i)
        if os.path.exists(training_completed_flag):
            logger.info(f"Training already completed for {training_path} and skipped.")
            return
        
        os.makedirs(training_path, exist_ok=True)
        self._merge_dataset(i)
        training_set_config = copy.deepcopy(self.pretrain_config["Datahub"])
        training_set_config["data_path"] = training_set_path
        validation_set_config = copy.deepcopy(self.pretrain_config["Datahub"])
        validation_set_config["data_path"] = validation_set_path
        new_trainer_config = {
            "Splitter": {
                "method": "random",
                "parts": [
                    {"name": "training", "dataset": "training"},
                    {"name": "validation", "dataset": "validation"}
                ],
                "save": False
            }
        }
        if self.continual_learning and i > 0:
            new_trainer_config.update({
                "resume": 2,
                "refresh_patience": True,
                "refresh_best_score": True
            })
        if self.reset_parameters and i <= 0:
            new_trainer_config["reset_parameters"] = True
        self._make_model_config(i + 1,
            dump_path=training_config_path, 
            new_config_path=self.training_config_path,
            new_datahub_config={
                "datasets": {
                    "training": training_set_config,
                    "validation": validation_set_config,
                },
                "global_transforms": self.pretrain_config["Datahub"]["transforms"]
            },
            new_trainer_config=new_trainer_config,
            model_params_updater=self._training_model_params_updater
        )
        enerzyme_subprocess = subprocess.Popen(
            ["enerzyme", "train", "-c", training_config_path, "-o", training_path],
            stdout=sys.stdout,
            stderr=sys.stderr
        )
        enerzyme_subprocess.wait()
        shutil.move(os.path.join(training_path, self._get_model_str(i + 1)), os.path.join(self.output_path))
        f = open(os.path.join(training_path, "training_completed"), "w")
        f.close()
        logger.info(f"Training finished for {training_path}")

    def launch(self) -> None:
        self._get_topology()
        self._get_model_config()
        initial_model_path = os.path.join(self.output_path, self._get_model_str(0))
        if self.pretrain_path is not None and not os.path.exists(initial_model_path):
            os.symlink(os.path.join(self.pretrain_path, self._get_model_str()), initial_model_path)

        self._init_structure_pool()

        for i in range(self.n_iterations):
            # simulation step
            self._prediction_step(i)
            self._simulation_step(i)
            self._collection_step(i)
            self._extraction_step(i)
            self._annotation_step(i)
            self._training_step(i)
