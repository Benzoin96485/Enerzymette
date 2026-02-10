from typing import List

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


def is_backbone(resname: str, atomname: str) -> bool:
    return resname in resname_list and atomname in ["N", "C", "CA", "O"]

def is_Calpha(resname: str, atomname: str) -> bool:
    return resname in resname_list and atomname == "CA"

def is_water_oxygen(resname: str, atomname: str) -> bool:
    return resname == "HOH" and atomname == "O"

def is_ion(resname: str, atomname: str) -> bool:
    return resname in ["MG", "NA", "CL", "K", "CA"]

index_type_map = {
    "backbone": is_backbone,
    "C_alpha": is_Calpha,
    "O_water": is_water_oxygen,
    "ions": is_ion
}

def get_indices(pdb_path: str, idx_start_from: int, index_types: List[str]) -> List[int]:
    indices = []
    atom_count = idx_start_from
    with open(pdb_path, "r") as f:
        for line in f.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                resname = line[17:20].strip()
                atomname = line[11:16].strip()
                for index_type in index_types:
                    if index_type_map[index_type](resname, atomname):
                        indices.append(atom_count)
                atom_count += 1
    return indices
