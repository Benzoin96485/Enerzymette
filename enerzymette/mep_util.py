from typing import List, Tuple
import re

def find_intermediate_indices(energies: List[float]) -> List[int]:
    if len(energies) < 3:
        return []
    intermediate_indices = []
    for i in range(1, len(energies) - 1):
        if energies[i] < energies[i-1] and energies[i] < energies[i+1]:
            intermediate_indices.append(i)
    return intermediate_indices

def find_ci_index(energies: List[float]) -> int:
    ci_index = -1
    max_energy = float("-inf")
    for i, energy in enumerate(energies):
        if energy > max_energy:
            max_energy = energy
            ci_index = i
    return ci_index


def _find_rate_determining_step(intermediate_indices, energies, ci_index) -> Tuple[int, int, int]:
    new_reactant_index = 0
    new_product_index = len(energies) - 1
    for intermediate_index in sorted(intermediate_indices):
        if intermediate_index < ci_index:
            new_reactant_index = intermediate_index
        elif intermediate_index > ci_index:
            new_product_index = intermediate_index
            break
    return new_reactant_index, new_product_index


def find_rate_determining_step(intermediate_indices, energies, ci_index) -> Tuple[int, int, int]:
    if ci_index < 0:
        ci_index = find_ci_index(energies)
    new_reactant_index, new_product_index = _find_rate_determining_step(intermediate_indices, energies, ci_index)
    if intermediate_indices and new_reactant_index == 0 and new_product_index == len(energies) - 1:
        ci_index = find_ci_index(energies)
        new_reactant_index, new_product_index = _find_rate_determining_step(intermediate_indices, energies, ci_index)
    if intermediate_indices and new_reactant_index == 0 and new_product_index == len(energies) - 1:
        ci_index = -1
    return new_reactant_index, new_product_index, ci_index

def letters_to_int(s: str) -> int:
    """
    Converts a string of letters to an integer, treating it as a 26-adic number.
    'a' is 1, 'b' is 2, ..., 'z' is 26, 'aa' is 27, etc.

    Args:
        s (str): The input string consisting of lowercase letters.

    Returns:
        int: The corresponding integer value.
    """
    result = 0
    for char in s:
        result = result * 26 + (ord(char) - ord('a') + 1)
    return result

def int_to_letters(n: int) -> str:
    """
    Converts an integer to a string of lowercase letters, treating it as a 26-adic number.
    1 -> 'a', 2 -> 'b', ..., 26 -> 'z', 27 -> 'aa', etc.

    Args:
        n (int): The integer to convert (should be >= 1).

    Returns:
        str: The corresponding string of lowercase letters.
    """
    if n < 1:
        raise ValueError("Input must be a positive integer")
    result = ''
    while n > 0:
        n -= 1
        result = chr(ord('a') + (n % 26)) + result
        n //= 26
    return result

def find_new_name(existing_names: List[str]) -> str:
    species_indices = []
    conformer_indices = []
    for name in existing_names:
        species_index, conformer_index = re.match(r"(\d+)([a-z]+)", name).groups()
        species_indices.append(int(species_index))
        conformer_indices.append(letters_to_int(conformer_index))
    max_species_index = max(species_indices)
    max_conformer_index = max(conformer_indices)
    new_species_index = max_species_index
    new_conformer_index = max_conformer_index + 1
    return f"{new_species_index}{int_to_letters(new_conformer_index)}"