import os
from typing import Optional, List

from ..logger import logger
from ..plumed_config_generator import get_plumed_patch
from .workflow import (
    build_enerzyme_simulate_cmd,
    copy_reaction_local_minima,
    find_lowest_local_minima,
    run_elementary_reaction_scan,
    run_scan_chain,
    write_standalone_scan_config,
)


class EnerzymeScanLauncher:
    def __init__(self,
        reactant_path: str,
        output_path: str,
        model_path: str,
        reference_path: str,
        model_config_path: Optional[str]=None,
        reactant_name: str="1a",
        product_name: str="2a",
        n_steps: int=25,
        plumed_patch_key: Optional[str]=None,
        plumed_cv_config: Optional[dict]=None,
    ):
        self.reactant_path = reactant_path
        self.output_path = output_path
        self.local_minima_path = os.path.join(self.output_path, "local_minima")
        self.model_path = model_path
        self.reference_path = reference_path
        self.model_config_path = model_config_path
        logger.info(f"Using model config: {self.model_config_path}")
        self.reactant_name = reactant_name
        self.product_name = product_name
        self.n_steps = n_steps
        self.plumed_patch_key = plumed_patch_key
        self.plumed_cv_config = plumed_cv_config or {}
        self.plumed_patch = get_plumed_patch(plumed_patch_key) if plumed_patch_key is not None else None
        from .io import infer_reference_type
        self.reference_type = infer_reference_type(reference_path)
        self.reference = self.parse_reference(reference_path, self.reference_type)
        self.constraint_freeze_xyz = self.reference.get("constraint_freeze", {}).get("xyz", [])
        self.constraint_scan = self.reference.get("constraint_scan", {})
        if plumed_patch_key is not None:
            logger.info(f"Using PLUMED CV-plugin scan mode (patch: {plumed_patch_key})")
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
        if reference_type == "scan_config":
            logger.info(f"Parsing system information from scan config {reference_path}")
            from .io import parse_scan_config
            if not os.path.exists(reference_path):
                logger.warning(f"Reference file {reference_path} does not exist")
                return {"main": {}, "constraint_freeze": {}, "constraint_scan": {}}
            return parse_scan_config(reference_path, self.output_path)
        raise NotImplementedError(f"Reference type {reference_type} is not supported")

    def copy_local_minima(self, reactant_name: str, product_name: str):
        elementary_reaction_path = os.path.join(
            self.output_path, f"{reactant_name}-{product_name}"
        )
        copy_reaction_local_minima(
            self.local_minima_path,
            reactant_name,
            product_name,
            elementary_reaction_path,
        )

    def find_lowest_local_minima(self):
        return find_lowest_local_minima(self.local_minima_path)

    def launch(self):
        run_scan_chain(
            self.output_path,
            self.reactant_path,
            run_elementary_reaction=self.launch_elementary_reaction,
            local_minima_path=self.local_minima_path,
            reactant_name=self.reactant_name,
            product_name=self.product_name,
            write_results=True,
            handle_errors=True,
        )

    def launch_elementary_reaction(self,
        reactant_path: str,
        elementary_reaction_path: str,
        target_value: Optional[float]=None,
        target_structure_path: Optional[str]=None,
    ):
        scan_task = "plumed_scan" if self.plumed_patch_key is not None else "scan"
        model_config_arg = ["-mc", self.model_config_path] if self.model_config_path is not None else []

        def write_config(task, initial_structure_path, config_path, **kwargs):
            self.write_config(
                task=task,
                initial_structure_path=initial_structure_path,
                config_path=config_path,
                target_value=kwargs.get("target_value"),
                target_structure_path=kwargs.get("target_structure_path"),
            )

        def build_simulate_cmd(config_path, output_path, with_plumed_patch):
            return self._enerzyme_simulate_cmd(
                config_path,
                output_path,
                model_config_arg,
                with_plumed_patch=with_plumed_patch,
            )

        return run_elementary_reaction_scan(
            reactant_path,
            elementary_reaction_path,
            write_config=write_config,
            build_simulate_cmd=build_simulate_cmd,
            scan_task=scan_task,
            target_value=target_value,
            target_structure_path=target_structure_path,
        )

    def _enerzyme_simulate_cmd(
        self,
        config_path: str,
        output_path: str,
        model_config_arg: List[str],
        with_plumed_patch: bool = False,
    ) -> List[str]:
        return build_enerzyme_simulate_cmd(
            config_path,
            output_path,
            self.model_path,
            plumed_patch=self.plumed_patch if with_plumed_patch else None,
            model_config_arg=model_config_arg,
        )

    def write_config(self,
        task: str,
        initial_structure_path: str,
        config_path: str,
        target_value: Optional[float]=None,
        target_structure_path: Optional[str]=None,
    ):
        idx_start_from = 0 if self.reference_type == "scan_config" else 1
        write_standalone_scan_config(
            config_path,
            task=task,
            initial_structure_path=initial_structure_path,
            charge=self.charge,
            multiplicity=self.multiplicity,
            constraint_freeze_xyz=self.constraint_freeze_xyz,
            idx_start_from=idx_start_from,
            constraint_scan=self.constraint_scan,
            plumed_patch_key=self.plumed_patch_key,
            plumed_cv_config=self.plumed_cv_config,
            n_steps=self.n_steps,
            target_value=target_value,
            target_structure_path=target_structure_path,
        )
