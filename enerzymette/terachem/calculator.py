import re
from pathlib import Path
from .ase_io import write_terachem, read_terachem_outputs, write_xyz, clean_scr
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

    def __init__(self, label: str='terachem', keep_mo: bool=True):
        super().__init__('terachem',
                         implemented_properties=['energy', 'forces', 'dipole'])

        self._label = label
        self.keep_mo = keep_mo
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

        new_scr_path = Path(parameters["scrdir"]) if "scrdir" in parameters else None
        if new_scr_path is not None:
            find_mo_flags = clean_scr(new_scr_path, keep_mo=self.keep_mo, label=self._label)
            for key, value in find_mo_flags.items():
                if value:
                    print(f"move molecular orbitals {key} from {new_scr_path} to {value}")
            coordinates_path, scr_path = write_terachem(directory / self.inputname, atoms, kw, find_mo_flags)
        else:
            coordinates_path, scr_path = write_terachem(directory / self.inputname, atoms, kw, False)
            clean_scr(scr_path, keep_mo=self.keep_mo, label=self._label)
        write_xyz(directory / coordinates_path, atoms)

    def read_results(self, directory):
        return read_terachem_outputs(directory, directory / self.outputname)

    def load_profile(self, cfg, **kwargs):
        return TerachemProfile.from_config(cfg, self.name, **kwargs)


class Terachem(GenericFileIOCalculator):
    """Class for doing TERACHEM calculations.

    Example:

      calc = Terachem(charge=0, mult=1, terachemblocks='''run gradient
basis 6-31gs
method b3lyp
end
''')
    """

    def __init__(self, *, profile=None, directory='.', label='terachem', keep_mo=True, **kwargs):
        """Construct Terachem-calculator object.

        Parameters
        ==========
        charge: int

        mult: int

        terachemblocks: str


        Examples
        ========
        Use default values:

        >>> from enerzymette.terachem.calculator import Terachem
        >>> He = Atoms(
        ...     'He',
        ...     calculator=Terachem(
        ...         charge=0,
        ...         mult=1,
        ...         terachemblocks='''run gradient
basis 6-31gs
method b3lyp
end
'''))
        """

        super().__init__(template=TerachemTemplate(label=label, keep_mo=keep_mo),
                         profile=profile, directory=directory,
                         parameters=kwargs)
