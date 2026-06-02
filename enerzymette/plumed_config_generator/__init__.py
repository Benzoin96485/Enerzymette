"""PLUMED collective-variable plugins for Enerzymette launchers and Enerzyme.

See ``README.md`` in this package for the plugin contract, YAML schema
(``task: plumed_scan``, ``sampling.cv: plumed``), and active-learning
``--initial-scan`` options.
"""

from __future__ import annotations

from dataclasses import dataclass
from importlib import import_module
from types import ModuleType
from typing import Callable, Dict, List, Optional, Tuple

import ase.io

from ._engine import (
    ReactionCoordinate,
    generate_scan_restraint,
    generate_steered_md,
)

__all__ = [
    "PLUMED_CV_PLUGINS",
    "PlumedCvPluginSpec",
    "ReactionCoordinate",
    "generate_scan_restraint",
    "generate_steered_md",
    "get_plumed_patch",
    "get_reaction_coordinate_fn",
    "get_scan_config_fn",
    "get_scan_config_generator_name",
    "get_steered_config_fn",
    "get_steered_config_generator_name",
    "list_plumed_cv_plugin_keys",
    "register_plumed_cv_plugin",
    "resolve_scan_endpoints",
]


@dataclass(frozen=True)
class PlumedCvPluginSpec:
    """Metadata for a built-in PLUMED CV plugin."""

    key: str
    module_name: str
    description: str = ""

    @property
    def reaction_coordinate_fn_name(self) -> str:
        return f"get_{self.key}_reaction_coordinate"

    @property
    def steered_config_fn_name(self) -> str:
        return f"get_{self.key}_config"

    @property
    def scan_config_fn_name(self) -> str:
        return f"get_{self.key}_scan_config"


# Registry of CV plugins shipped with Enerzymette. Keys are used by CLI
# (-pp sammt) and by scantoolkit / altoolkit when emitting Enerzyme configs.
PLUMED_CV_PLUGINS: Dict[str, PlumedCvPluginSpec] = {
    "sammt": PlumedCvPluginSpec(
        key="sammt",
        module_name=".sammt",
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
    description: str = "",
    overwrite: bool = False,
) -> None:
    """Register a CV plugin at runtime (e.g. from an extension package).

    ``module_name`` is either a dotted import path or a relative name under
    ``enerzymette.plumed_config_generator`` (e.g. ``".my_cv"``).
    The target module must implement the function names implied by ``key``;
    see ``README.md``.
    """
    if key in PLUMED_CV_PLUGINS and not overwrite:
        raise ValueError(f"PLUMED CV plugin {key!r} is already registered")
    PLUMED_CV_PLUGINS[key] = PlumedCvPluginSpec(
        key=key, module_name=module_name, description=description
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


def _get_plugin_callable(key: str, fn_name: str) -> Callable:
    module = _load_plugin_module(key)
    fn = getattr(module, fn_name, None)
    if fn is None:
        raise AttributeError(
            f"PLUMED CV plugin {key!r} ({module.__file__}) has no function {fn_name!r}"
        )
    return fn


def get_plumed_patch(key: str) -> str:
    """Path to the plugin module file for Enerzyme ``simulate -pp``."""
    return _load_plugin_module(key).__file__


def get_reaction_coordinate_fn(key: str) -> Callable:
    return _get_plugin_callable(key, _plugin_spec(key).reaction_coordinate_fn_name)


def get_steered_config_fn(key: str) -> Callable:
    return _get_plugin_callable(key, _plugin_spec(key).steered_config_fn_name)


def get_scan_config_fn(key: str) -> Callable:
    return _get_plugin_callable(key, _plugin_spec(key).scan_config_fn_name)


def get_steered_config_generator_name(key: str) -> str:
    """``plumed_config_generator.name`` for steered MD (``task: plumed``)."""
    return _plugin_spec(key).steered_config_fn_name


def get_scan_config_generator_name(key: str) -> str:
    """``plumed_config_generator.name`` for flexible scan (``task: plumed_scan``)."""
    return _plugin_spec(key).scan_config_fn_name


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
    rc_fn = get_reaction_coordinate_fn(plumed_patch_key)
    rc = rc_fn(system, idx_start_from=idx_start_from, **plumed_cv_config)
    x0 = rc.initial_value
    if target_value is not None:
        x1 = target_value
    elif target_structure_path is not None:
        target_system = ase.io.read(target_structure_path, index=-1)
        target_rc = rc_fn(target_system, idx_start_from=idx_start_from, **plumed_cv_config)
        x1 = target_rc.initial_value
    else:
        dist_to_lower = abs(x0 - rc.lower_bound)
        dist_to_upper = abs(x0 - rc.upper_bound)
        x1 = rc.lower_bound if dist_to_lower >= dist_to_upper else rc.upper_bound
    return x0, x1, num, rc
