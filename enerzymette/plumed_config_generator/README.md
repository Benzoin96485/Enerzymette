# PLUMED Config Generators

Enerzymette builds PLUMED input through class-based generator plugins. A plugin is a Python module that exposes a `PlumedConfigGenerator` subclass. Enerzyme receives the plugin module through `enerzyme simulate -pp <path>`, instantiates the configured class with the current `ase.Atoms` object plus YAML parameters, and calls the configured method to produce `plumed.dat`.

The design separates three layers:

- `PlumedConfigGenerator` in `_engine.py`: system-independent workflow methods for steered MD, naive steered MD, scan restraints, and PLUMED unit preamble.
- Chemistry-specific generator subclasses such as `SAMMTConfigGenerator` in `sammt.py`: atom discovery and main reaction-coordinate definition.
- Proton-transfer plugins in `proton_transfer.py`, such as `local_opes`: optional auxiliary CV/bias builders that can be inserted into any generator through `proton_transfer`.

## User Interface

A simulation YAML selects a generator class and one method:

```yaml
Simulation:
  task: plumed
  idx_start_from: 1
  plumed_config_generator:
    name: SAMMTConfigGenerator
    method: standard_steered_md
  integrate:
    integrator: Langevin
    n_step: 100000
    time_step: 0.5
    temperature_in_K: 500
  sampling:
    params:
      plumed_config:
        dump_interval: 20
        lower_bound: -2
        upper_bound: 2
        reference_pdb_file: cluster-capped.pdb
        substrate: G
        nucleophile: "O2'"
System:
  structure_file: initial.xyz
```

Enerzyme passes these values to `SAMMTConfigGenerator(system, ...)`, then calls `standard_steered_md(...)`. Supported base methods are:

- `standard_steered_md`: starts from the current main reaction coordinate and runs one round trip across `[lower_bound, upper_bound]`.
- `naive_steered_md`: first pulls the coordinate to the nearest bound over `warmup_steps`, then pulls to the opposite bound.
- `scan`: applies a static PLUMED `RESTRAINT` at `target_value`; Enerzyme supplies `target_value` for each `task: plumed_scan` point. Scan configs do not insert proton-transfer plugins, even if `proton_transfer` is present in `plumed_config`.

For scan jobs, use the same generator class with `method: scan`:

```yaml
Simulation:
  task: plumed_scan
  idx_start_from: 1
  optimize:
    optimizer: LBFGS
  plumed_config_generator:
    name: SAMMTConfigGenerator
    method: scan
  sampling:
    cv: plumed
    params:
      x0: 0.4
      x1: -1.4
      num: 25
      plumed_config:
        dump_interval: 20
        lower_bound: -2
        upper_bound: 2
        reference_pdb_file: cluster-capped.pdb
        substrate: G
        nucleophile: "O2'"
```

## SAMMT Generator

`SAMMTConfigGenerator` is for SAM-dependent methyltransferase systems. Its main coordinate is:

```text
dd = d(CE, nucleophile) - d(SD, CE)
```

Indices can come from a reference PDB:

```yaml
plumed_config:
  reference_pdb_file: cluster-capped.pdb
  substrate: G
  nucleophile: "O2'"
```

Or from explicit atom indices in the same convention as `Simulation.idx_start_from`:

```yaml
plumed_config:
  index_sulphur: 292
  index_methyl_carbon: 293
  index_nucleophile: 388
```

The descriptive index names exposed by `get_indices()` are `sulphur`, `sulfur`, `methyl_carbon`, and `nucleophile`. Proton-transfer configs can refer to these names, for example `donor: nucleophile`.

## Proton Transfer

Proton-transfer support lives outside `_engine.py` in `proton_transfer.py`. The base generator only asks that module whether a configured plugin should append extra PLUMED lines.

Enable proton-transfer OPES by adding a nested `proton_transfer` mapping under `plumed_config`:

```yaml
plumed_config:
  dump_interval: 20
  lower_bound: -2
  upper_bound: 2
  reference_pdb_file: cluster-capped.pdb
  substrate: G
  nucleophile: "O2'"
  proton_transfer:
    enabled: true
    plugin: local_opes
    donor: nucleophile
    flavor: nearest_distance
    scope_file: structure_pool/000.pt_scope.json
    state_file: opes_state.data
    restart: false
    topology_mol_file: cluster.mol
    opes_barrier: 20
```

`local_opes` resolves a donor, donor-bound transfer proton(s), and nearby N/O acceptors. If `scope_file` exists, it reuses the saved scope; otherwise it resolves a new scope and writes it. If required atoms cannot be found, the generator warns and emits the main-coordinate PLUMED config without proton-transfer bias.

Important `proton_transfer` fields:

- `enabled`: switch for inserting proton-transfer lines.
- `plugin`: proton-transfer plugin key; currently `local_opes`.
- `donor`: descriptive atom name from `get_indices()` or an explicit atom index. SAMMT defaults this to `nucleophile` when omitted.
- `proton`: optional explicit proton selector or selector list.
- `acceptor`: optional explicit acceptor selector or selector list.
- `flavor`: `nearest_distance`, `coordination`, `grouped_acceptor`, or `coupled_dd_pt`. The default is `nearest_distance`.
- `scope_file`: optional JSON sidecar for resolved donor/proton/acceptor indices.
- `restart`: add PLUMED `RESTART` and read `state_file` as `STATE_RFILE`.
- `state_file`: OPES state file written by `STATE_WFILE`.
- `topology_mol_file`: optional RDKit-readable mol file used to filter chemically saturated N/O acceptors.
- `opes_pace`, `opes_barrier`, `opes_biasfactor`, `opes_state_wstride`: OPES tuning.
- `geometry_walls`, `max_donor_h_distance`, `max_donor_acceptor_distance`: optional walls used by nontrivial flavors.

The default `nearest_distance` flavor defines a distance-gap CV from the selected donor-bound proton to the nearest acceptor. `coordination` uses coordination numbers, `grouped_acceptor` compares acceptor groups, and `coupled_dd_pt` combines proton transfer with the main `dd` coordinate and donor-acceptor contact.

OPES requires a PLUMED build with the OPES module enabled, for example by loading an OPES-enabled PLUMED module or setting `PLUMED_KERNEL`.

## Developer Interface

To add a generator, create a module with a `PlumedConfigGenerator` subclass:

```python
from enerzymette.plumed_config_generator import PlumedConfigGenerator

class MyConfigGenerator(PlumedConfigGenerator):
    default_print_args = "rc,mr.*"

    def get_indices(self):
        return {"atom_a": self.index_a, "atom_b": self.index_b}

    def define_main_rc(self):
        return "rc", "rc: DISTANCE ATOMS=1,2 NOPBC"

    def calc_main_rc(self):
        return self.system.get_distance(0, 1, mic=False)
```

The base class constructor already records:

- `system`: the current runtime `ase.Atoms` object.
- `idx_start_from`: index convention used by user configs.
- `preamble`: PLUMED unit preamble, defaulting to Angstrom/fs/kJ-mol units.
- `reference_pdb`: path from `reference_pdb` or `reference_pdb_file`.
- `reference_mol`: an RDKit mol object, loaded from `reference_mol_file` or `topology_mol_file` when available.
- `proton_transfer_config`: parsed `ProtonTransferConfig`.

Subclass constructors should accept chemistry-specific keyword arguments, call `super().__init__(system, **kwargs)`, resolve any atom indices, and store them on `self`.

Register built-in or external generators with the registry in `__init__.py`:

```python
register_plumed_cv_plugin(
    "mykey",
    "my_package.my_plumed_module",
    class_name="MyConfigGenerator",
)
```

Launchers use the registry to write:

```yaml
plumed_config_generator:
  name: MyConfigGenerator
  method: standard_steered_md
```

The plugin module path for Enerzyme is resolved with `get_plumed_patch("mykey")`.

## Proton-Transfer Plugin Interface

To add another proton-transfer enhanced-sampling strategy, implement `ProtonTransferPlugin` and register an instance:

```python
from enerzymette.plumed_config_generator import (
    ProtonTransferConfig,
    ProtonTransferPlugin,
    register_proton_transfer_plugin,
)

class MyProtonTransferPlugin(ProtonTransferPlugin):
    name = "my_pt"

    def append_to_plumed(
        self,
        generator,
        plumed_config,
        config: ProtonTransferConfig,
        dump_interval: int,
        integrate_config: dict,
    ):
        indices = generator.get_indices()
        lines = list(plumed_config)
        # Insert custom PT CV/bias lines before PRINT and update PRINT args.
        return lines

register_proton_transfer_plugin(MyProtonTransferPlugin())
```

The plugin receives the fully initialized `PlumedConfigGenerator`, so it can use `generator.system`, `generator.idx_start_from`, `generator.reference_pdb`, `generator.reference_mol`, `generator.get_indices()`, and the main reaction coordinate definition if needed. Users select it with:

```yaml
proton_transfer:
  enabled: true
  plugin: my_pt
```

`proton_transfer.py` also exposes reusable helpers used by `local_opes`, including `ProtonTransferScope`, `append_proton_transfer_to_plumed()`, `generate_proton_transfer_cv_lines()`, and `generate_opes_explore()`.

## Launcher Integration

`resolve_scan_endpoints()` instantiates the registered generator and uses `build_reaction_coordinate()` to compute the current coordinate. It chooses scan endpoints from an explicit target value, a target structure, or the farther configured bound.

`altoolkit` and `scantoolkit` emit class-based scan configs automatically when a PLUMED plugin key is supplied. Active learning also injects per-structure `proton_transfer.scope_file`, `proton_transfer.state_file`, `proton_transfer.restart`, and `proton_transfer.topology_mol_file` when proton transfer is enabled in the base YAML.

## Debugging

Installed packages may be non-editable. During development, prepend source trees to `PYTHONPATH` or reinstall editable packages:

```bash
export PYTHONPATH=/path/to/Enerzymette:/path/to/Enerzyme:$PYTHONPATH
```

Use the `enerzyme-dev` conda environment for smoke tests in this workspace.
