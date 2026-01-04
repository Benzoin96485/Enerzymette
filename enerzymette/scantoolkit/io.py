from ..terachem.io import parse_terachem_input, write_terachem_input
import ase.io

def update_terachem_scan_input(terachem_input_file: str, updated_structure_xyz: str, output_path: str) -> None:
    info = parse_terachem_input(terachem_input_file)
    i0 = info["constraint_scan"]["bond"]["i0"]
    i1 = info["constraint_scan"]["bond"]["i1"]
    updated_structure = ase.io.read(updated_structure_xyz, index=-1)
    x0 = updated_structure.get_distance(i0 - 1, i1 - 1)
    info["constraint_scan"]["bond"]["x0"] = x0
    write_terachem_input(info, output_path)