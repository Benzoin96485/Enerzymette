import re
from .io import write_terachem, read_terachem_outputs, write_xyz, clean_scr
from ase.calculators.genericfileio import (
    BaseProfile,
    CalculatorTemplate,
    GenericFileIOCalculator,
)


def get_version_from_terachem_header(terachem_header):
    match = re.search(r'TeraChem version (\S+)', terachem_header, re.M)
    return match.group(1)


class TerachemProfile(BaseProfile):
    def version(self):
        # XXX Allow MPI in argv; the version call should not be parallel.
        from ase.calculators.genericfileio import read_stdout
        stdout = read_stdout([self.command, "-v"])
        return get_version_from_terachem_header(stdout)

    def get_calculator_command(self, inputfile):
        return [inputfile]


class TerachemTemplate(CalculatorTemplate):
    _label = 'terachem'

    def __init__(self):
        super().__init__('terachem',
                         implemented_properties=['energy', 'forces', 'dipole'])

        self.inputname = f'{self._label}.inp'
        self.outputname = f'{self._label}.out'
        self.errorname = f'{self._label}.err'

    def execute(self, directory, profile) -> None:
        profile.run(directory, self.inputname, self.outputname,
                    errorfile=self.errorname)

    def write_input(self, profile, directory, atoms, parameters, properties):
        parameters = dict(parameters)

        kw = dict()
        kw.update(parameters)

        coordinates_path, scr_path = write_terachem(directory / self.inputname, atoms, kw)
        clean_scr(scr_path)
        write_xyz(directory / coordinates_path, atoms)

    def read_results(self, directory):
        return read_terachem_outputs(directory, directory / self.outputname)

    def load_profile(self, cfg, **kwargs):
        return TerachemProfile.from_config(cfg, self.name, **kwargs)


class Terachem(GenericFileIOCalculator):
    """Class for doing TERACHEM calculations.

    Example:

      calc = Terachem(charge=0, mult=1, orcasimpleinput='B3LYP def2-TZVP',
        orcablocks='%pal nprocs 16 end')
    """

    def __init__(self, *, profile=None, directory='.', **kwargs):
        """Construct ORCA-calculator object.

        Parameters
        ==========
        charge: int

        mult: int

        orcasimpleinput : str

        orcablocks: str


        Examples
        ========
        Use default values:

        >>> from ase.calculators.orca import ORCA
        >>> h = Atoms(
        ...     'H',
        ...     calculator=ORCA(
        ...         charge=0,
        ...         mult=1,
        ...         directory='water',
        ...         orcasimpleinput='B3LYP def2-TZVP',
        ...         orcablocks='%pal nprocs 16 end'))

        """

        super().__init__(template=TerachemTemplate(),
                         profile=profile, directory=directory,
                         parameters=kwargs)
