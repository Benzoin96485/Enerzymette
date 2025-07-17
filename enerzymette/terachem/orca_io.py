from typing import Dict, Any, List, Optional
from .io import parse_terachem_input, write_terachem_input, run_terachem, read_energy, read_grad, clean_scr


def parse_orca_extinp(orca_extinp_file: str) -> Dict[str, Any]:
    if not orca_extinp_file.endswith(".extinp.tmp"):
        raise ValueError(f"Unsupported format: {orca_extinp_file}")
    basename = orca_extinp_file.replace(".extinp.tmp", "")
    with open(orca_extinp_file, "r") as f:
        lines = f.readlines()
    xyz_file = lines[0].split("#")[0].strip()
    charge = int(lines[1].split("#")[0].strip())
    multiplicity = int(lines[2].split("#")[0].strip())
    return {
        "xyz_file": xyz_file,
        "charge": charge,
        "multiplicity": multiplicity,
        "basename": basename,
    }


def write_orca_output(results: Dict[str, Any], basename: str):
    with open(f"{basename}.engrad", "w") as f:
        f.write(f'''
#
# Number of atoms: must match the XYZ
#
{len(results["outputs"]["Fa"])}
#
# The current total energy in Eh
#
{results["outputs"]["E"] / results["units"]["Hartree_in_E"]}
#
# The current gradient in Eh/bohr: Atom1X, Atom1Y, Atom1Z, Atom2X, etc.
#
''')
        for atom_Fa in results["outputs"]["Fa"]:
            for component in atom_Fa:
                f.write(f"{-component / results['units']['Hartree_in_E'] * results['units']['Bohr_in_R']}\n")


def write_orca_engrad(basename: str, energy: float, forces: List[List[float]], N: Optional[int] = None):
    with open(f"{basename}.engrad", "w") as f:
        f.write(f'''
#
# Number of atoms: must match the XYZ
#
{N if N is not None else len(forces)}
#
# The current total energy in Eh
#
{energy:.12f}
#
# The current gradient in Eh/bohr: Atom1X, Atom1Y, Atom1Z, Atom2X, etc.
#
''')
        for atom_force in forces:
            for component in atom_force:
                f.write(f"{component:.12f}\n")


def run(orca_extinp_file: str, terachem_input_template: str):
    orca_extinp_info = parse_orca_extinp(orca_extinp_file)
    basename = orca_extinp_info["basename"]
    scr_path = "./scr_" + basename
    find_mo_flags = clean_scr(scr_path=scr_path, keep_mo=True, basename=basename)
    terachem_input_info = parse_terachem_input(terachem_input_template)
    main_info = terachem_input_info["main"]
    main_info["run"] = "gradient"
    main_info["coordinates"] = orca_extinp_info["xyz_file"]
    main_info["charge"] = orca_extinp_info["charge"]
    main_info["spinmult"] = orca_extinp_info["multiplicity"]
    main_info.pop("guess", None)
    if find_mo_flags["c0"]:
        main_info["guess"] = find_mo_flags["c0"]
    elif find_mo_flags["ca"] and find_mo_flags["cb"]:
        main_info["guess"] = (find_mo_flags["ca"], find_mo_flags["cb"])
    main_info["scrdir"] = scr_path
    write_terachem_input(info={"main": main_info}, terachem_input_file=basename + ".terachem.inp")
    output_file, error_file = run_terachem(terachem_input_file=basename + ".terachem.inp", basename=basename)
    with open(output_file, "r") as f:
        output_lines = f.readlines()
    energy = read_energy(output_lines)
    forces = read_grad(output_lines)
    write_orca_engrad(basename, energy=energy, forces=forces)
