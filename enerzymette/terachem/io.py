from typing import List, Optional, Dict, Any
from collections import defaultdict
from shutil import rmtree
import subprocess
import os


def read_energy(lines: List[str]) -> Optional[float]:
    """Read energy."""
    energy = None
    for line in lines:
        if 'FINAL ENERGY:' in line:
            try:
                energy = float(line.split()[-2])
            except:
                energy = float('nan')
    return energy
    

def read_grad(lines: List[str]) -> Optional[List[List[float]]]:
    """Read gradients."""
    gradients = None
    start_read_grad = False
    read_grad_idx = 0
    for line in lines:
        if 'Gradient units are Hartree/Bohr' in line:
            start_read_grad = True
            gradients = []
            continue
        if read_grad_idx > 2 and "--" in line:
            start_read_grad = False
            read_grad_idx = 0
            continue
        if start_read_grad:
            read_grad_idx += 1
            if read_grad_idx > 2:
                gx, gy, gz = line.split()
                gradients.append([float(gx), float(gy), float(gz)])
    return gradients


def parse_terachem_input(terachem_input_file: str) -> Dict[str, Dict[str, Any]]:
    with open(terachem_input_file, "r") as f:
        lines = f.readlines()

    # parse the main section
    info = {"main": {}}
    current_section = "main"
    for line in lines:
        clean_line = line.split("#")[0].strip()
        if not clean_line:
            continue
        if clean_line == "end":
            if current_section == "main":
                current_section = None
            else:
                raise ValueError(f"Unexpected end in section: {current_section}")
        elif current_section is None and clean_line.startswith("$"):
            current_section = clean_line[:-1]
            info[current_section] = defaultdict(list)
        elif clean_line == "$end":
            if current_section == "main" or current_section is None:
                raise ValueError(f"Unexpected $end in section: {current_section}")
            else:
                current_section = None
        elif current_section == "main":
            key, *value = clean_line.split()
            if len(value) == 1:
                info["main"][key] = value[0]
            else:
                info[current_section][key] = value
        elif current_section is not None:
            info[current_section]["content"].append(clean_line)

    return info


def write_terachem_input(info: Dict[str, Any], terachem_input_file: str):
    with open(terachem_input_file, "w") as f:
        main_section = info["main"]
        for key, value in main_section.items():
            if isinstance(value, tuple):
                f.write(f"{key} {' '.join(value)}\n")
            else:
                f.write(f"{key} {value}\n")
        f.write("end\n")
        for section, content in info.items():
            if section == "main":
                continue
            f.write(f"${section}\n")
            for line in content["content"]:
                f.write(f"{line}\n")
            f.write("$end\n")


def clean_scr(scr_path: str, keep_mo: bool, basename: Optional[str]=None):
    if basename is None:
        basename = os.path.basename(scr_path).lstrip("scr").lstrip("_")
    # for some unknown reason, terachem will create a scr directory in lower case
    real_scr_path = scr_path.lower()
    find_mo_flags = {
        "c0": '',
        "ca": '',
        "cb": '',
    }
    if os.path.exists(real_scr_path):
        if keep_mo:
            for mo_name in ["c0", "ca", "cb"]:
                mo_path = os.path.join(real_scr_path, mo_name)
                if os.path.exists(mo_path):
                    new_mo_path = f"{mo_name}_{basename}"
                    find_mo_flags[mo_name] = new_mo_path
                    os.rename(mo_path, new_mo_path)
        rmtree(real_scr_path)
    return find_mo_flags


def run_terachem(terachem_input_file: str, basename: Optional[str]=None):
    if basename is None:
        result = subprocess.run(["terachem", terachem_input_file], capture_output=True, text=True)
        return result.stdout, result.stderr
    else:
        output_file = f"{basename}.out"
        error_file = f"{basename}.err"
        with open(output_file, "w") as out_fd, open(error_file, "w") as err_fd:
            result = subprocess.run(["terachem", terachem_input_file], text=True, stdout=out_fd, stderr=err_fd)
        return output_file, error_file
    
