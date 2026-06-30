"""PLUMED collective-variable plugins for Enerzymette launchers and Enerzyme.

See ``README.md`` in this package for the plugin contract, YAML schema
(``task: plumed_scan``, ``sampling.cv: plumed``), and active-learning
``--initial-scan`` options.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from types import ModuleType
from typing import Dict, List, Optional, Tuple, Type

import ase.io

from ._engine import (
    PlumedConfigGenerator,
    ReactionCoordinate,
    generate_scan_restraint,
    generate_steered_md,
)
from .proton_transfer import (
    LocalOpesProtonTransferPlugin,
    PROTON_TRANSFER_PLUGINS,
    ProtonTransferConfig,
    ProtonTransferPlugin,
    ProtonTransferScope,
    append_optional_proton_transfer,
    build_proton_transfer_config,
    register_proton_transfer_plugin,
)

__all__ = [
    "PLUMED_CV_PLUGINS",
    "PlumedConfigGenerator",
    "PlumedCvPluginSpec",
    "LocalOpesProtonTransferPlugin",
    "PROTON_TRANSFER_PLUGINS",
    "ProtonTransferConfig",
    "ProtonTransferPlugin",
    "ProtonTransferScope",
    "append_optional_proton_transfer",
    "build_proton_transfer_config",
    "ReactionCoordinate",
    "generate_scan_restraint",
    "generate_steered_md",
    "get_config_generator_class",
    "get_config_generator_name",
    "get_plumed_patch",
    "get_scan_method_name",
    "get_steered_method_name",
    "list_plumed_cv_plugin_keys",
    "register_plumed_cv_plugin",
    "register_proton_transfer_plugin",
    "resolve_scan_endpoints",
]


@dataclass(frozen=True)
class PlumedCvPluginSpec:
    """Metadata for a built-in PLUMED CV plugin."""

    key: str
    module_name: str
    class_name: str
    steered_method_name: str = "standard_steered_md"
    scan_method_name: str = "scan"
    description: str = ""


# Registry of CV plugins shipped with Enerzymette. Keys are used by CLI
# (-pp sammt) and by scantoolkit / altoolkit when emitting Enerzyme configs.
PLUMED_CV_PLUGINS: Dict[str, PlumedCvPluginSpec] = {
    "sammt": PlumedCvPluginSpec(
        key="sammt",
        module_name=".sammt",
        class_name="SAMMTConfigGenerator",
        description=(
            "SAM methyltransferase coordinate dd = d(S–CE) − d(CE–Nu); "
            "indices from PDB (SAM SD/CE + substrate nucleophile) or explicit indices."
        ),
    ),
}


def register_plumed_cv_plugin(
    key: str,
    module_name: str,
    *,
    class_name: str,
    steered_method_name: str = "standard_steered_md",
    scan_method_name: str = "scan",
    description: str = "",
    overwrite: bool = False,
) -> None:
    """Register a CV plugin at runtime (e.g. from an extension package).

    ``module_name`` is either a dotted import path or a relative name under
    ``enerzymette.plumed_config_generator`` (e.g. ``".my_cv"``).
    The target module must expose ``class_name`` as a ``PlumedConfigGenerator``
    subclass; see ``README.md``.
    """
    if key in PLUMED_CV_PLUGINS and not overwrite:
        raise ValueError(f"PLUMED CV plugin {key!r} is already registered")
    PLUMED_CV_PLUGINS[key] = PlumedCvPluginSpec(
        key=key,
        module_name=module_name,
        class_name=class_name,
        steered_method_name=steered_method_name,
        scan_method_name=scan_method_name,
        description=description,
    )


def list_plumed_cv_plugin_keys() -> List[str]:
    return sorted(PLUMED_CV_PLUGINS.keys())


def _plugin_spec(key: str) -> PlumedCvPluginSpec:
    try:
        return PLUMED_CV_PLUGINS[key]
    except KeyError:
        known = ", ".join(list_plumed_cv_plugin_keys()) or "(none)"
        raise ValueError(f"Unknown PLUMED CV plugin key: {key!r}. Known keys: {known}") from None


def _load_plugin_module(key: str) -> ModuleType:
    spec = _plugin_spec(key)
    if spec.module_name.startswith("."):
        return import_module(spec.module_name, package=__name__)
    return import_module(spec.module_name)


def get_config_generator_class(key: str) -> Type[PlumedConfigGenerator]:
    spec = _plugin_spec(key)
    module = _load_plugin_module(key)
    cls = getattr(module, spec.class_name, None)
    if cls is None:
        raise AttributeError(
            f"PLUMED CV plugin {key!r} ({module.__file__}) has no class "
            f"{spec.class_name!r}"
        )
    if not issubclass(cls, PlumedConfigGenerator):
        raise TypeError(
            f"PLUMED CV plugin {key!r} class {spec.class_name!r} must inherit "
            "PlumedConfigGenerator"
        )
    return cls


def get_plumed_patch(key: str) -> str:
    """Path to the plugin module file for Enerzyme ``simulate -pp``."""
    return _load_plugin_module(key).__file__


def get_config_generator_name(key: str) -> str:
    return _plugin_spec(key).class_name


def get_steered_method_name(key: str) -> str:
    return _plugin_spec(key).steered_method_name


def get_scan_method_name(key: str) -> str:
    return _plugin_spec(key).scan_method_name


def resolve_scan_endpoints(
    system,
    idx_start_from: int,
    plumed_patch_key: str,
    plumed_cv_config: dict,
    num: int,
    target_value: Optional[float] = None,
    target_structure_path: Optional[str] = None,
) -> Tuple[float, float, int, ReactionCoordinate]:
    """Pick scan endpoints x0/x1 from a structure and :class:`ReactionCoordinate` bounds."""
    generator_cls = get_config_generator_class(plumed_patch_key)
    generator = generator_cls(system, idx_start_from=idx_start_from, **plumed_cv_config)
    rc = generator.build_reaction_coordinate(**plumed_cv_config)
    x0 = rc.initial_value
    if target_value is not None:
        x1 = target_value
    elif target_structure_path is not None:
        target_system = ase.io.read(target_structure_path, index=-1)
        target_generator = generator_cls(
            target_system,
            idx_start_from=idx_start_from,
            **plumed_cv_config,
        )
        target_rc = target_generator.build_reaction_coordinate(**plumed_cv_config)
        x1 = target_rc.initial_value
    else:
        dist_to_lower = abs(x0 - rc.lower_bound)
        dist_to_upper = abs(x0 - rc.upper_bound)
        x1 = rc.lower_bound if dist_to_lower >= dist_to_upper else rc.upper_bound
    return x0, x1, num, rc
