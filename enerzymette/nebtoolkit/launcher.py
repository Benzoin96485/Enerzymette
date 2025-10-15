from typing import Optional, Literal, Dict, Any
from threading import Thread
import ase.io
from ase.units import kcal, mol, Ha
import os, subprocess, time, sys, traceback, shutil
from glob import glob
from .logger import logger
from .io import write_orca_neb_in, make_backup, read_energy, redirect_output, parse_neb_csv, get_mep_path_info
from .analysis import find_intermediate_indices, find_rate_determining_step, find_new_name, check_neb_convergence


class EnerzymeNEBLauncher:
    def __init__(self, 
        reactant_path: str,
        product_path: str,
        output_path: str,
        model_path: str,
        reference_path: str,
        server_config_path: str,
        n_images: int=25,
        port: int=5000,
        reactant_name: str="1a",
        product_name: str="2a",
        interrupt_strategy: Literal["none", "stdout"]="stdout"
    ):
        # sanity check
        self.orca_exe = os.environ.get("ORCA_PATH", None)
        if self.orca_exe is None:
            raise RuntimeError("ORCA_PATH is not set")
        self.n_images = n_images
        self.reactant_path = reactant_path
        self.product_path = product_path
        self.output_path = output_path
        self.local_minima_path = os.path.join(self.output_path, "local_minima")
        self.model_path = model_path
        self.server_config_path = server_config_path
        self.port = port
        self.reactant_name = reactant_name
        self.product_name = product_name
        self.reference = self.parse_reference(reference_path, "terachem_input")
        self.constraint_freeze_xyz = self.reference.get("constraint_freeze", {}).get("xyz", [])
        logger.info(f"Constraint freeze xyz: {self.constraint_freeze_xyz}")
        self.charge = self.reference.get("main", {}).get("charge", 0)
        logger.info(f"Charge: {self.charge}")
        self.multiplicity = self.reference.get("main", {}).get("spinmult", 1)
        logger.info(f"Multiplicity: {self.multiplicity}")
        self.interrupt_strategy = interrupt_strategy
        self.enerzyme_subprocess = None

    def copy_local_minima(self, reactant_name: str, product_name: str):
        elementary_reaction_path = os.path.join(self.output_path, f"{reactant_name}-{product_name}")
        reactant_path = os.path.join(elementary_reaction_path, f"neb_reactant.xyz")
        product_path = os.path.join(elementary_reaction_path, f"neb_product.xyz")
        if os.path.exists(reactant_path):
            shutil.copy(reactant_path, os.path.join(self.local_minima_path, "reactant", f"{reactant_name}.xyz"))
        else:
            logger.warning(f"Reactant file does not exist for finished reaction {reactant_name}-{product_name}")
        if os.path.exists(product_path):
            shutil.copy(product_path, os.path.join(self.local_minima_path, "product", f"{product_name}.xyz"))
        else:
            logger.warning(f"Product file does not exist for finished reaction {reactant_name}-{product_name}")

    def find_lowest_local_minima(self):
        reactant_paths = glob(os.path.join(self.local_minima_path, "reactant", "*.xyz"))
        product_paths = glob(os.path.join(self.local_minima_path, "product", "*.xyz"))
        reactant_energies = [read_energy(reactant_path) for reactant_path in reactant_paths]
        product_energies = [read_energy(product_path) for product_path in product_paths]
        lowest_reactant_energy = min(reactant_energies)
        lowest_product_energy = min(product_energies)
        lowest_reactant_name = reactant_paths[reactant_energies.index(lowest_reactant_energy)]
        lowest_product_name = product_paths[product_energies.index(lowest_product_energy)]
        return lowest_reactant_name, lowest_product_name, lowest_reactant_energy, lowest_product_energy

    def launch(self):
        # launch enerzyme server
        enerzyme_subprocess = subprocess.Popen(
            ["enerzyme", "listen", "-c", self.server_config_path, "-o", self.output_path, "-m", self.model_path, "-b", f"0.0.0.0:{self.port}"], 
            cwd=self.output_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT
        )
        self.enerzyme_subprocess = enerzyme_subprocess

        # monitor stdout, wait until the server is ready ("waitress | Serving" displayed in the file)
        while True:
            line = self.enerzyme_subprocess.stdout.readline()
            if line and "waitress | Serving" in line.decode("utf-8"):
                break
        logger.info(f"Enerzyme server is ready")

        os.makedirs(self.local_minima_path, exist_ok=True)
        os.makedirs(os.path.join(self.local_minima_path, "reactant"), exist_ok=True)
        os.makedirs(os.path.join(self.local_minima_path, "product"), exist_ok=True)

        try:
            reactant_names = [self.reactant_name]
            product_names = [self.product_name]
            current_reactant_name = self.reactant_name
            current_product_name = self.product_name
            current_reactant_path = self.reactant_path
            current_product_path = self.product_path
            current_ts_path = None

            resume = False
            existing_reactions = []
            csv_path = os.path.join(self.output_path, "neb.csv")
            if os.path.exists(csv_path):
                existing_reactions = parse_neb_csv(csv_path)
            if len(existing_reactions) > 0:
                for reactant_name, product_name in existing_reactions[:-1]:
                    self.copy_local_minima(reactant_name, product_name)

                if len(existing_reactions) > 1:
                    for reactant_name, product_name in existing_reactions[1:]:
                        reactant_names.append(reactant_name)
                        product_names.append(product_name)
                current_reactant_name, current_product_name = existing_reactions[-1]
                logger.info(f"Restarting from last reaction {current_reactant_name}-{current_product_name}")
                csv_fp = open(csv_path, "a")
                resume = True
            else:
                csv_fp = open(csv_path, "w")

            while True:
                reaction_name = f"{current_reactant_name}-{current_product_name}"
                elementary_reaction_path = os.path.join(self.output_path, reaction_name)
                if not resume:
                    logger.info(f"Creating reaction {reaction_name}")

                elementary_reaction_info = self.check_elementary_reaction(elementary_reaction_path)
                if elementary_reaction_info is None:
                    if resume:
                        logger.info(f"Rerun reaction {reaction_name}")
                        current_reactant_path = os.path.join(elementary_reaction_path, f"reactant.xyz")
                        current_product_path = os.path.join(elementary_reaction_path, f"product.xyz")
                        current_ts_path = os.path.join(elementary_reaction_path, f"ts.xyz")
                        if not os.path.exists(current_reactant_path):
                            current_ts_path = None
                    else:
                        csv_fp.write(f"{current_reactant_name},{current_product_name}\n")
                    elementary_reaction_info = self.launch_elementary_reaction(current_reactant_path, current_product_path, elementary_reaction_path, current_ts_path)
                else:
                    logger.info(f"Information collected for last reaction {reaction_name}")

                self.copy_local_minima(current_reactant_name, current_product_name)
                resume = False

                if elementary_reaction_info["intermediate_indices"]:
                    mep_path_info = elementary_reaction_info["mep_path_info"]
                    reactant_index, product_index, ts_index = find_rate_determining_step(elementary_reaction_info["intermediate_indices"], mep_path_info["energies"])
                    current_reactant_path = os.path.join(elementary_reaction_path, f"{reaction_name}-{reactant_index}.xyz")
                    current_product_path = os.path.join(elementary_reaction_path, f"{reaction_name}-{product_index}.xyz")
                    current_ts_path = os.path.join(elementary_reaction_path, f"{reaction_name}-{ts_index}.xyz")
                    with open(current_reactant_path, "w") as f:
                        f.writelines(mep_path_info["xyzblocks"][reactant_index])
                    with open(current_product_path, "w") as f:
                        f.writelines(mep_path_info["xyzblocks"][product_index])
                    with open(current_ts_path, "w") as f:
                        f.writelines(mep_path_info["xyzblocks"][ts_index])
                    if reactant_index != 0:
                        current_reactant_name = find_new_name(reactant_names)
                        reactant_names.append(current_reactant_name)
                        logger.info(f"New reactant {current_reactant_name} of rate determining step is found from image {reactant_index} of reaction {reaction_name}")
                    if product_index != self.n_images - 1:
                        current_product_name = find_new_name(product_names)
                        product_names.append(current_product_name)
                        logger.info(f"New product {current_product_name} of rate determining step is found from image {product_index} of reaction {reaction_name}")
                    logger.info(f"TS initial guess for new reaction {current_reactant_name}-{current_product_name} is found from image {ts_index} of reaction {reaction_name}")
                elif not elementary_reaction_info["converged"]:
                    logger.warning(f"NEB did not converge for reaction {reaction_name}, restarting...")
                    self.restart_elementary_reaction(elementary_reaction_path)
                else:
                    logger.info(f"NEB converged for reaction {reaction_name}!")
                    ts_energy = read_energy(os.path.join(elementary_reaction_path, "neb_NEB-CI_converged.xyz"))
                    lowest_reactant_name, lowest_product_name, lowest_reactant_energy, lowest_product_energy = self.find_lowest_local_minima()
                    energy_span = (ts_energy - lowest_reactant_energy) * Ha / (kcal / mol)
                    energy_change = (lowest_product_energy - lowest_reactant_energy) * Ha / (kcal / mol)
                    logger.info(f"Reaction energy span: {energy_span:.2f} kcal/mol between reactant {lowest_reactant_name} and transition state of reaction {reaction_name}")
                    logger.info(f"Reaction energy change: {energy_change:.2f} kcal/mol between reactant {lowest_reactant_name} and product {lowest_product_name}")
                    break
            
            # clean up
            logger.info(f"Calculation finished. Terminating enerzyme server...")
        except Exception as e:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback_strs = traceback.format_exception(exc_type, exc_value, exc_traceback)
            logger.error(f"Error: {exc_type.__name__}: {exc_value}")
            logger.error(f"Error traceback: {"".join(traceback_strs)}")
            logger.error(f"Terminating enerzyme server...")

        csv_fp.close()
        subprocess.Popen(["enerzyme", "kill", "-u", f"0.0.0.0:{self.port}"], cwd=self.output_path)
        time.sleep(10)
        if self.enerzyme_subprocess.poll() is None:
            self.enerzyme_subprocess.kill()
        logger.info(f"Enerzyme server terminated")
        
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

    def write_enerzyme_wrapper(self, wrapper_path: str):
        wrapper_str = f"""enerzyme request -i $1 -u 0.0.0.0:{self.port}
"""
        with open(wrapper_path, "w") as f:
            f.write(wrapper_str)
        # set executable permission
        os.system(f"chmod +x {wrapper_path}")

    def check_elementary_reaction(self, elementary_reaction_path: str):
        mep_path_info = get_mep_path_info(elementary_reaction_path)
        if mep_path_info is not None:
            intermediate_indices = find_intermediate_indices(mep_path_info["energies"])
            neb_converged = check_neb_convergence(elementary_reaction_path)
            return {
                "mep_path_info": mep_path_info,
                "intermediate_indices": intermediate_indices,
                "converged": neb_converged
            }
        else:
            return None

    def monitor_elementary_reaction(self, elementary_reaction_path: str, neb_in_path: str):
        logger.info(f"Launching ORCA to run NEB with input file {neb_in_path}")
        neb_out_fp = open(os.path.join(elementary_reaction_path, "neb.out"), "w")
        neb_err_fp = open(os.path.join(elementary_reaction_path, "neb.err"), "w")
        enerzyme_out_fp = open("enerzyme.out", "w")
        enerzyme_out_thread = Thread(target=redirect_output, args=(self.enerzyme_subprocess.stdout, enerzyme_out_fp), daemon=True)
        enerzyme_out_thread.start()
        intermediate_indices = []
        orca_subprocess = subprocess.Popen(
            [self.orca_exe, "neb.in"], cwd=elementary_reaction_path,
            stdout=subprocess.PIPE,
            stderr=neb_err_fp
        )
        try:
            converged = False
            # neb_out_monitor_fp = open(os.path.join(elementary_reaction_path, "neb.out"), "r")
            if self.interrupt_strategy == "none":
                # wait for orca subprocess to finish
                orca_subprocess.wait()
                mep_path_info = get_mep_path_info(elementary_reaction_path)
                intermediate_indices = find_intermediate_indices(mep_path_info["energies"])
            elif self.interrupt_strategy == "stdout":
        # Read ORCA subprocess stdout line by line and kill upon detecting intermediate minimum
                while orca_subprocess.poll() is None:
                    line = orca_subprocess.stdout.readline()
                    if line:
                        line_decoded = line.decode("utf-8")
                        if "Possible intermediate minimum found at" in line_decoded:
                            logger.info(f"ORCA warned: {line_decoded.strip()}, checking if there is an intermediate minimum...")
                            mep_path_info = get_mep_path_info(elementary_reaction_path)
                            intermediate_indices = find_intermediate_indices(mep_path_info["energies"])
                            logger.info(f"Intermediate indices in the chain: {intermediate_indices}")
                            if len(intermediate_indices) > 0:
                                neb_out_fp.flush()
                                orca_subprocess.terminate()
                                break
                        elif "THE NEB OPTIMIZATION HAS CONVERGED" in line_decoded:
                            converged = True
                        neb_out_fp.write(line_decoded)
                        neb_out_fp.flush()

        except Exception as e:
            neb_out_fp.close()
            neb_err_fp.close()
            enerzyme_out_fp.close()
            enerzyme_out_thread.join()
            raise e
        return {
            "converged": converged,
            "intermediate_indices": intermediate_indices,
            "mep_path_info": mep_path_info
        }

    def launch_elementary_reaction(self, 
        reactant_path: str, 
        product_path: str, 
        elementary_reaction_path: str,
        ts_path: Optional[str]=None,
    ) -> Dict[str, Any]:
        # create a folder
        if os.path.exists(elementary_reaction_path):
            logger.warning(f"Output path {elementary_reaction_path} already exists, it could be overwritten")
        else:
            os.makedirs(elementary_reaction_path)
        
        # read reactant and product
        logger.info(f"Reading reactant from {reactant_path}")
        reactant_atoms = ase.io.read(reactant_path, index=-1)
        ase.io.write(os.path.join(elementary_reaction_path, f"reactant.xyz"), reactant_atoms, format="xyz")
        logger.info(f"Reading product from {product_path}")
        product_atoms = ase.io.read(product_path, index=-1)
        ase.io.write(os.path.join(elementary_reaction_path, f"product.xyz"), product_atoms, format="xyz")
        logger.info(f"Reactant and product written to {elementary_reaction_path}")

        # read ts
        if ts_path is not None:
            logger.info(f"Reading guessed transition state from {ts_path}")
            if os.path.exists(ts_path):
                ts_atoms = ase.io.read(ts_path, index=-1)
                ase.io.write(os.path.join(elementary_reaction_path, f"ts.xyz"), ts_atoms, format="xyz")
                logger.info(f"Guessed transition state written to {elementary_reaction_path}")
            else:
                logger.warning(f"Guessed transition state file {ts_path} does not exist")

        # write enerzyme wrapper
        wrapper_path = os.path.join(elementary_reaction_path, "enerzyme_wrapper.sh")
        self.write_enerzyme_wrapper(wrapper_path)

        # write neb.in
        neb_in_path = os.path.join(elementary_reaction_path, "neb.in")
        write_orca_neb_in(
            neb_in_path, wrapper_path, 
            n_images=self.n_images,
            constraint_freeze_xyz=self.constraint_freeze_xyz, 
            charge=self.charge, 
            multiplicity=self.multiplicity
        )
        
        # spawn orca subprocess to run neb
        return self.monitor_elementary_reaction(elementary_reaction_path, neb_in_path)

    def restart_elementary_reaction(self, 
        elementary_reaction_path: str
    ):
        neb_in_path = os.path.join(elementary_reaction_path, "neb.in")
        wrapper_path = os.path.join(elementary_reaction_path, "enerzyme_wrapper.sh")
        for filename in ["neb_MEP.allxyz", "neb.out", "neb.err"]:
            make_backup(elementary_reaction_path, filename)
        
        write_orca_neb_in(
            neb_in_path, wrapper_path, 
            n_images=self.n_images,
            restart=True,
            pre_opt=False,
            use_ts=False,
            constraint_freeze_xyz=self.constraint_freeze_xyz, 
            charge=self.charge, 
            multiplicity=self.multiplicity,
        )
        return self.monitor_elementary_reaction(elementary_reaction_path, neb_in_path)
