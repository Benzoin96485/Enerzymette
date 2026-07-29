import os, subprocess, pickle, copy, sys, shutil, glob, json
from random import shuffle
from typing import Tuple, Optional, Dict, Callable, Literal, List
import yaml
import ase.io
from ..external_calculator import get_calculator_patch
from ..plumed_config_generator import get_plumed_patch, get_steered_method_name
from ..scantoolkit.workflow import (
    build_enerzyme_simulate_cmd,
    run_elementary_reaction_scan,
    run_scan_chain,
    write_base_scan_task_config,
)
from ..logger import logger
from .get_index import get_indices
from .structures_manifest import SystemManifestEntry, load_systems_manifest


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
        initial_structures_config_path: Optional[str]=None,
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
        self.initial_structures_config_path = initial_structures_config_path
        self.multi_system_mode = initial_structures_config_path is not None
        self.systems_manifest: Optional[List[SystemManifestEntry]] = None
        self._validate_launch_options()

    def _validate_launch_options(self) -> None:
        if not self.multi_system_mode:
            return

        if self.initial_scan:
            raise NotImplementedError(
                "Multi-system active learning does not support --initial-scan"
            )
        if self.reference_pdb_path is not None:
            raise ValueError(
                "Use reference_pdb in --initial-structures-config instead of "
                "-rp/--reference_pdb_path"
            )
        if self.template_sdf_path is not None:
            raise ValueError(
                "Use reference_sdf in --initial-structures-config instead of "
                "-ts/--template_sdf_path"
            )
        if self.initial_xyz_path is not None:
            raise ValueError(
                "Use reference_xyz in --initial-structures-config instead of "
                "-ix/--initial_xyz_path"
            )

        self.systems_manifest = load_systems_manifest(self.initial_structures_config_path)
        logger.info(
            f"Loaded multi-system manifest with {len(self.systems_manifest)} systems"
        )

    def _topology_dir_for_system(self, name: str) -> str:
        return os.path.join(self.output_path, "topology", name)

    def _topology_paths_for_system(self, name: str) -> Tuple[str, str, str]:
        topo_dir = self._topology_dir_for_system(name)
        return (
            os.path.join(topo_dir, "cluster.png"),
            os.path.join(topo_dir, "cluster.xyz"),
            os.path.join(topo_dir, "cluster.mol"),
        )

    def _cluster_mol_path_for_system(self, name: str) -> str:
        return self._topology_paths_for_system(name)[2]

    def _run_enerzyme_bond(
        self,
        reference_pdb: str,
        output_img_path: str,
        output_mol_path: str,
        reference_sdf: Optional[str] = None,
    ) -> None:
        enerzyme_bond_args = [
            "enerzyme", "bond",
            "-p", reference_pdb,
            "-i", output_img_path,
            "-m", output_mol_path,
        ]
        if reference_sdf is not None:
            enerzyme_bond_args.extend(["-t", reference_sdf])
        enerzyme_subprocess = subprocess.Popen(enerzyme_bond_args)
        enerzyme_subprocess.wait()

    def _apply_simulation_topology(
        self,
        simulation_config: Dict,
        *,
        initial_xyz_path: Optional[str] = None,
        charge: Optional[int] = None,
        reference_pdb_path: Optional[str] = None,
        backbone_indices: Optional[List[int]] = None,
        Calpha_indices: Optional[List[int]] = None,
    ) -> Dict:
        if initial_xyz_path is not None:
            simulation_config["System"]["structure_file"] = initial_xyz_path
        if charge is not None:
            simulation_config["System"]["charge"] = charge
        if backbone_indices is not None and Calpha_indices is not None:
            idx_start_from = simulation_config["Simulation"]["idx_start_from"]
            if self.restraint_mode == "hard":
                simulation_config["Simulation"]["constraint"] = {
                    "fix_atom": {
                        "indices": [idx + idx_start_from for idx in backbone_indices]
                    }
                }
            elif self.restraint_mode == "soft":
                Hookean_k = simulation_config["Simulation"]["constraint"].get(
                    "Hookean_allpairs", {}
                ).get("k", 0.05)
                simulation_config["Simulation"]["constraint"] = {
                    "Hookean_allpairs": {
                        "indices": [idx + idx_start_from for idx in Calpha_indices],
                        "k": Hookean_k,
                    }
                }
        if reference_pdb_path is not None:
            simulation_config["Simulation"]["sampling"]["params"]["plumed_config"][
                "reference_pdb_file"
            ] = reference_pdb_path
        return simulation_config

    def _get_multi_system_topology(self) -> None:
        from rdkit.Chem import MolFromMolFile, MolToXYZFile, GetFormalCharge

        for system in self.systems_manifest:
            topo_dir = self._topology_dir_for_system(system.name)
            os.makedirs(topo_dir, exist_ok=True)
            img_path, xyz_path, mol_path = self._topology_paths_for_system(system.name)
            if os.path.exists(mol_path):
                logger.info(
                    f"Topology for {system.name} already exists at {topo_dir}, skipping bond."
                )
            else:
                self._run_enerzyme_bond(
                    system.reference_pdb,
                    img_path,
                    mol_path,
                    reference_sdf=system.reference_sdf,
                )
                logger.info(
                    f"Generated topology for {system.name} under {topo_dir} "
                    f"(cluster.mol, cluster.png)"
                )

            cluster_mol = MolFromMolFile(mol_path, removeHs=False)
            charge = GetFormalCharge(cluster_mol)
            initial_xyz_path = system.reference_xyz
            if not system.reference_xyz_explicit:
                MolToXYZFile(cluster_mol, xyz_path)
                initial_xyz_path = xyz_path
                system.reference_xyz = xyz_path
                system.source_structure = xyz_path

            backbone_indices = get_indices(system.reference_pdb, 0, ["backbone"])
            Calpha_indices = get_indices(
                system.reference_pdb, 0, ["C_alpha", "O_water"]
            )
            system.charge = charge
            system.backbone_indices = backbone_indices
            system.Calpha_indices = Calpha_indices
            logger.info(
                f"Using total charge ({charge}) of the cluster parsed from the "
                f"reference structure: {system.reference_pdb} with the small molecule "
                f"template: {system.reference_sdf}"
            )
            if not system.reference_xyz_explicit:
                logger.info(
                    f"Using initial structure: {initial_xyz_path} converted from the "
                    f"reference structure: {system.reference_pdb}"
                )
            logger.info(
                f"Using topology information {mol_path} parsed from the reference "
                f"structure: {system.reference_pdb} and the small molecule template: "
                f"{system.reference_sdf}"
            )

    def _get_topology(self) -> None:
        if self.multi_system_mode:
            self._get_multi_system_topology()
            return
        if self.reference_pdb_path is not None:
            from rdkit.Chem import MolFromMolFile, MolToXYZFile, GetFormalCharge
            self._run_enerzyme_bond(
                self.reference_pdb_path,
                self.output_img_path,
                self.output_mol_path,
                reference_sdf=self.template_sdf_path,
            )
            cluster_mol = MolFromMolFile(self.output_mol_path, removeHs=False)
            self.charge = GetFormalCharge(cluster_mol)
            if self.initial_xyz_path is None:
                MolToXYZFile(cluster_mol, self.output_xyz_path)
                self.initial_xyz_path = self.output_xyz_path

            self.backbone_indices = get_indices(self.reference_pdb_path, 0, ["backbone"])
            self.Calpha_indices = get_indices(
                self.reference_pdb_path, 0, ["C_alpha", "O_water"]
            )
            logger.info(f"Using total charge ({self.charge}) of the cluster parsed from the reference structure: {self.reference_pdb_path} with the small molecule template: {self.template_sdf_path}")
            logger.info(f"Using initial structure: {self.initial_xyz_path} converted from the reference structure: {self.reference_pdb_path}")
            
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

    def _pt_scope_path_for_entry(self, entry_path: str) -> str:
        base, _ = os.path.splitext(entry_path)
        return f"{base}.pt_scope.json"

    def _pt_state_path_for_entry(self, entry_path: str) -> str:
        base, _ = os.path.splitext(entry_path)
        return f"{base}.opes_state.data"

    def _plumed_config_from_simulation(self, simulation_config: Dict) -> Optional[Dict]:
        return (
            simulation_config.get("Simulation", {})
            .get("sampling", {})
            .get("params", {})
            .get("plumed_config")
        )

    def _get_structure_pool_entry(self, iteration: int) -> Tuple[int, Dict]:
        entries = self.structure_pool_state["entries"]
        pool_idx = iteration % len(entries)
        return pool_idx, entries[pool_idx]

    def _simulation_config_path_for_pool_entry(self, entry: Dict) -> str:
        return entry.get("simulation_config", self.simulation_config_path)

    def _system_manifest_entry_for_pool_entry(
        self, entry: Dict
    ) -> Optional[SystemManifestEntry]:
        if not self.multi_system_mode or not entry.get("name"):
            return None
        return next(
            (
                system
                for system in self.systems_manifest
                if system.name == entry["name"]
            ),
            None,
        )

    def _apply_topology_for_pool_entry(
        self, simulation_config: Dict, entry: Dict
    ) -> Dict:
        if self.multi_system_mode:
            system = self._system_manifest_entry_for_pool_entry(entry)
            if system is not None:
                return self._apply_simulation_topology(
                    simulation_config,
                    initial_xyz_path=system.reference_xyz,
                    charge=system.charge,
                    reference_pdb_path=system.reference_pdb,
                    backbone_indices=system.backbone_indices,
                    Calpha_indices=system.Calpha_indices,
                )
            return simulation_config

        if self.reference_pdb_path is not None:
            return self._apply_simulation_topology(
                simulation_config,
                initial_xyz_path=self.initial_xyz_path,
                charge=self.charge,
                reference_pdb_path=self.reference_pdb_path,
                backbone_indices=self.backbone_indices,
                Calpha_indices=self.Calpha_indices,
            )
        return simulation_config

    def _load_base_simulation_config(self) -> Dict:
        with open(self.simulation_config_path, "r") as handle:
            simulation_config = yaml.load(handle, Loader=yaml.FullLoader)
        return self._apply_topology_for_pool_entry(simulation_config, {})

    def _load_simulation_config_for_entry(self, entry: Dict) -> Dict:
        config_path = self._simulation_config_path_for_pool_entry(entry)
        with open(config_path, "r") as handle:
            simulation_config = yaml.load(handle, Loader=yaml.FullLoader)
        return self._apply_topology_for_pool_entry(simulation_config, entry)

    def _inject_proton_transfer_config(
        self,
        iteration: int,
        simulation_config: Dict,
        simulation_path: str,
    ) -> Dict:
        plumed_config = self._plumed_config_from_simulation(simulation_config)
        if not plumed_config or not plumed_config.get("proton_transfer"):
            return simulation_config
        proton_transfer_config = plumed_config["proton_transfer"]
        if proton_transfer_config is True:
            proton_transfer_config = {"enabled": True}
            plumed_config["proton_transfer"] = proton_transfer_config

        pool_idx, entry = self._get_structure_pool_entry(iteration)
        scope_file = self._pt_scope_path_for_entry(entry["path"])
        state_pool_file = self._pt_state_path_for_entry(entry["path"])

        proton_transfer_config["scope_file"] = scope_file
        state_sim_file = os.path.join(simulation_path, "opes_state.data")
        proton_transfer_config["state_file"] = state_sim_file
        if os.path.exists(self.output_mol_path):
            proton_transfer_config["topology_mol_file"] = self.output_mol_path

        if entry.get("run") and os.path.exists(state_pool_file):
            proton_transfer_config["restart"] = True
            shutil.copy(state_pool_file, state_sim_file)
            logger.info(
                f"Iteration {iteration}: OPES restart for pool entry {pool_idx} "
                f"from {state_pool_file}"
            )
        else:
            proton_transfer_config["restart"] = False

        return simulation_config

    def _persist_proton_transfer_state(
        self,
        iteration: int,
        simulation_path: str,
        simulation_config: Dict,
    ) -> None:
        plumed_config = self._plumed_config_from_simulation(simulation_config)
        if not plumed_config or not plumed_config.get("proton_transfer"):
            return

        pool_idx, entry = self._get_structure_pool_entry(iteration)
        state_pool_file = self._pt_state_path_for_entry(entry["path"])
        plumed_config = self._plumed_config_from_simulation(simulation_config) or {}
        proton_transfer_config = plumed_config.get("proton_transfer")
        if proton_transfer_config is True:
            proton_transfer_config = {}
        state_sim_file = (proton_transfer_config or {}).get(
            "state_file", os.path.join(simulation_path, "opes_state.data")
        )
        if os.path.exists(state_sim_file):
            shutil.copy(state_sim_file, state_pool_file)
            logger.info(
                f"Iteration {iteration}: saved OPES state for pool entry {pool_idx} "
                f"to {state_pool_file}"
            )

    def _enerzyme_simulate_cmd(
        self,
        config_path: str,
        output_path: str,
        model_config_path: str,
        with_plumed_patch: bool = False,
    ) -> List[str]:
        return build_enerzyme_simulate_cmd(
            config_path,
            output_path,
            self.output_path,
            calculator_patch=self.calculator_patch,
            plumed_patch=self.plumed_patch if with_plumed_patch else None,
            model_config_path=model_config_path,
        )

    def _run_initial_scan(self) -> None:
        initial_scan_path = os.path.join(self.output_path, "initial_scan")
        local_minima_path = os.path.join(initial_scan_path, "local_minima")
        completed_flag = os.path.join(initial_scan_path, "initial_scan_completed")
        if os.path.exists(completed_flag):
            logger.info(f"Initial scan already completed at {initial_scan_path}, skipping.")
            return

        model_config_path = os.path.join(initial_scan_path, "model_config.yaml")
        if not os.path.exists(model_config_path):
            self._make_model_config(0, model_config_path)

        base_config = self._load_base_simulation_config()

        def write_config(task, structure_path, config_path, **kwargs):
            write_base_scan_task_config(
                config_path,
                task=task,
                structure_path=structure_path,
                base_config=base_config,
                plumed_patch_key=self.plumed_patch_key if task == "plumed_scan" else None,
                n_steps=self.n_initial_scan_steps,
                scan_target_value=kwargs.get("target_value"),
                scan_target_structure_path=kwargs.get("target_structure_path"),
                yaml_flow_style=False,
            )

        def build_simulate_cmd(config_path, output_path, with_plumed_patch):
            return build_enerzyme_simulate_cmd(
                config_path,
                output_path,
                self.output_path,
                calculator_patch=self.calculator_patch,
                plumed_patch=self.plumed_patch if with_plumed_patch else None,
                model_config_path=model_config_path,
            )

        def run_elementary(reactant_path, elementary_reaction_path, target_structure_path=None):
            return run_elementary_reaction_scan(
                reactant_path,
                elementary_reaction_path,
                write_config=write_config,
                build_simulate_cmd=build_simulate_cmd,
                scan_task="plumed_scan",
                target_structure_path=target_structure_path,
                redirect_stdio=True,
            )

        run_scan_chain(
            initial_scan_path,
            base_config["System"]["structure_file"],
            run_elementary_reaction=run_elementary,
            local_minima_path=local_minima_path,
            log_prefix="Initial scan: ",
        )

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
            if self.multi_system_mode:
                updated = False
                for entry in self.structure_pool_state["entries"]:
                    if "cluster_mol_path" not in entry and entry.get("name"):
                        entry["cluster_mol_path"] = self._cluster_mol_path_for_system(
                            entry["name"]
                        )
                        updated = True
                    system = self._system_manifest_entry_for_pool_entry(entry)
                    if system is not None:
                        for key, value in (
                            ("simulation_config", system.simulation_config),
                            ("reference_pdb", system.reference_pdb),
                            ("source_structure", system.source_structure),
                            ("charge", system.charge),
                            ("backbone_indices", system.backbone_indices),
                            ("Calpha_indices", system.Calpha_indices),
                        ):
                            if value is not None and entry.get(key) != value:
                                entry[key] = value
                                updated = True
                if updated:
                    self._save_structure_pool()
            logger.info(
                f"Loaded structure pool with {len(self.structure_pool_state['entries'])} entries"
            )
            return

        if self.initial_scan:
            self._run_initial_scan()
            self._build_structure_pool_from_initial_scan()
            return

        if self.multi_system_mode:
            pool_dir = self._structure_pool_dir()
            os.makedirs(pool_dir, exist_ok=True)
            pool_entries = []
            for pool_idx, system in enumerate(self.systems_manifest):
                dst = os.path.join(pool_dir, f"{pool_idx:03d}.xyz")
                shutil.copy(system.reference_xyz, dst)
                entry = {
                    "path": dst,
                    "run": False,
                    "name": system.name,
                    "simulation_config": system.simulation_config,
                    "reference_pdb": system.reference_pdb,
                    "source_structure": system.source_structure,
                    "cluster_mol_path": self._cluster_mol_path_for_system(system.name),
                    "charge": system.charge,
                    "backbone_indices": system.backbone_indices,
                    "Calpha_indices": system.Calpha_indices,
                }
                if system.reference_sdf:
                    entry["reference_sdf"] = system.reference_sdf
                pool_entries.append(entry)
            self.structure_pool_state = {"entries": pool_entries}
            self._save_structure_pool()
            logger.info(
                f"Structure pool initialized with {len(pool_entries)} systems "
                "from manifest"
            )
            return

        simulation_config = self._load_base_simulation_config()
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
        pool_idx, entry = self._get_structure_pool_entry(i)

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
                presimulation_config["Simulation"]["integrate"]["n_step"] = (
                    self.n_presimulation_steps_per_iteration
                )
                presimulation_config["System"]["structure_file"] = entry["path"]
                generator_cfg = (
                    presimulation_config["Simulation"].get("plumed_config_generator")
                    or {}
                )
                steered_method = get_steered_method_name(self.plumed_patch_key)
                use_restrained_plumed = generator_cfg.get("method") in {
                    steered_method,
                    "standard_steered_md",
                }
                if use_restrained_plumed:
                    # Equilibrate with additional restraints (e.g. UPPER_WALLS)
                    # but without MOVINGRESTRAINT.
                    presimulation_config["Simulation"]["task"] = "plumed"
                    generator_cfg = dict(generator_cfg)
                    generator_cfg["method"] = "standard_restrained_md"
                    presimulation_config["Simulation"][
                        "plumed_config_generator"
                    ] = generator_cfg
                    traj_path = os.path.join(simulation_path, "plumed.traj.xyz")
                    with_plumed_patch = True
                else:
                    presimulation_config["Simulation"].pop("sampling", None)
                    presimulation_config["Simulation"]["task"] = "md"
                    traj_path = presimulation_trajectory_path
                    with_plumed_patch = False
                with open(presimulation_config_path, "w") as f:
                    yaml.dump(presimulation_config, f, default_flow_style=False)
                subprocess.Popen(
                    self._enerzyme_simulate_cmd(
                        presimulation_config_path,
                        simulation_path,
                        simulation_model_config_path,
                        with_plumed_patch=with_plumed_patch,
                    ),
                    stdout=sys.stdout,
                    stderr=sys.stderr,
                ).wait()
                if os.path.exists(traj_path):
                    open(presimulation_completed_flag, "w").close()
                    # Restrained plumed writes plumed.traj.xyz; move aside so
                    # the subsequent steered MD can write a fresh trajectory.
                    if with_plumed_patch and traj_path != presimulation_trajectory_path:
                        if os.path.exists(presimulation_trajectory_path):
                            os.remove(presimulation_trajectory_path)
                        os.rename(traj_path, presimulation_trajectory_path)
            last_frame = ase.io.read(presimulation_trajectory_path, index=-1)
            ase.io.write(entry["path"], last_frame, format="extxyz")
            ase.io.write(initial_structure_path, last_frame, format="extxyz")
            self._save_structure_pool()
            logger.info(
                f"Iteration {i}: presimulation from pool entry {pool_idx}, "
                f"updated entry at {entry['path']}"
            )
            return

        if self.multi_system_mode:
            shutil.copy(entry["path"], initial_structure_path)
            logger.info(
                f"Iteration {i}: reusing pool entry {pool_idx} structure "
                f"from {entry['path']}"
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

        pool_idx, entry = self._get_structure_pool_entry(i)
        simulation_config = self._load_simulation_config_for_entry(entry)
        if i == 0:
            simulation_config["Simulation"].pop("uncertainty_calculator", None)
        else:
            average_E_var = self._get_max_E_var(i)
            if self.multi_system_mode:
                na = len(ase.io.read(entry["path"]))
                entry_name = entry.get("name", str(pool_idx))
                logger.info(
                    f"Pool entry {pool_idx} ({entry_name}) has {na} atoms."
                )
            else:
                if self.NA is None:
                    initial_structure = ase.io.read(
                        simulation_config["System"]["structure_file"]
                    )
                    self.NA = len(initial_structure)
                    logger.info(f"Initial structure has {self.NA} atoms.")
                na = self.NA
            simulation_config["Simulation"]["uncertainty_calculator"]["params"]["B"] = (
                average_E_var / na
            )

        simulation_config["System"]["structure_file"] = initial_structure_path
        simulation_config = self._inject_proton_transfer_config(
            i, simulation_config, simulation_path
        )

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
            # Intentionally do not write steered-MD end frames back into the
            # structure pool. Pool geometries are advanced only by
            # pre-simulation in _prepare_iteration_initial_structure. With
            # n_presimulation_steps_per_iteration == 0, reused multi-system
            # entries keep their seed/equilibrated coordinates.
        self._persist_proton_transfer_state(i, simulation_path, simulation_config)
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
        if self.multi_system_mode:
            _, entry = self._get_structure_pool_entry(i)
            cluster_mol_path = entry.get("cluster_mol_path")
            if cluster_mol_path:
                extraction_config["Extractor"]["reference_mol_path"] = cluster_mol_path
        elif self.reference_pdb_path is not None:
            extraction_config["Extractor"]["reference_mol_path"] = self.output_mol_path
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
