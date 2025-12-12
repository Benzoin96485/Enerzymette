from typing import List, Tuple
import re, os

def check_neb_convergence(elementary_reaction_path: str) -> bool:
    output_path = os.path.join(elementary_reaction_path, "neb.out")
    with open(output_path, "r") as f:
        for line in f:
            if "THE NEB OPTIMIZATION HAS CONVERGED" in line:
                return True
    return False

def check_latest_ci_index(elementary_reaction_path: str, ci_neb_pattern: re.Pattern) -> int:
    ci_index = -1
    output_path = os.path.join(elementary_reaction_path, "neb.out")
    with open(output_path, "r") as f:
        for line in f:
            match = ci_neb_pattern.match(line)
            if match is not None:
                ci_index = int(match.groups()[0])
    return ci_index
