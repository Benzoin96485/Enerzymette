from typing import Optional
from ase import Atoms
import ase.io
from ase.mep.neb import NEB
from ase.constraints import FixAtoms


def parse_constraints(constraints_file: str) -> list[FixAtoms]:
    '''
    Parse constraints from a terachem input file.
    
    Expected format:
    $constraint_freeze
    xyz 1,2,3,4          # GLN8
    xyz 19,20,21,22      # TYR9
    $end
    
    Returns a list of FixAtoms constraints with 0-based indices.
    '''
    all_indices = []
    
    with open(constraints_file, 'r') as f:
        lines = f.readlines()
    
    in_constraint_section = False
    
    for line in lines:
        line = line.strip()
        
        # Check for start of constraint section
        if line == '$constraint_freeze':
            in_constraint_section = True
            continue
        
        # Check for end of constraint section
        if line == '$end':
            in_constraint_section = False
            continue
        
        # Parse constraint lines within the section
        if in_constraint_section and line.startswith('xyz'):
            # Remove 'xyz' prefix and any trailing comments
            parts = line[3:].strip()
            if '#' in parts:
                parts = parts.split('#')[0].strip()
            
            # Parse comma-separated indices
            if parts:
                indices = [int(idx.strip()) for idx in parts.split(',')]
                # Convert to 0-based indexing (Terachem uses 1-based)
                indices_0based = [idx - 1 for idx in indices]
                all_indices.extend(indices_0based)
    
    # Remove duplicates and sort
    unique_indices = sorted(list(set(all_indices)))
    
    # Return a list of FixAtoms constraints
    return [FixAtoms(indices=unique_indices)]


def idpp(reactant_path: str, product_path: str, output_path: str, n_images: int, constraints_file: Optional[str]=None) -> None:
    '''
    IDPP interpolation between two molecules.
    '''
    # sanity check
    reactant: Atoms = ase.io.read(reactant_path, index=-1)
    product: Atoms = ase.io.read(product_path, index=-1)
    assert n_images > 2, "n_images must be greater than 2"

    images = [reactant]
    for i in range(n_images - 2):
        images.append(reactant.copy())
    images.append(product)

    if constraints_file is not None:
        constraints = parse_constraints(constraints_file)
        for image in images:
            image.set_constraint(constraints)

    neb = NEB(images)
    neb.interpolate(
        method='idpp',
        apply_constraint=True
    )
    ase.io.write(output_path, images, append=False)