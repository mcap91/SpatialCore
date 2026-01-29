"""
SpatialCore: Standardized spatial statistics tools for computational biology.

A thin, robust wrapper around standard libraries to ensure Python and R users
get the exact same result for the same biological question.
"""

__version__ = "0.4.5"

# Track which modules are available in this installation
_available_modules: list[str] = []


def _try_import(module_name: str):
    """Attempt to import a module, returning None if unavailable."""
    try:
        import importlib

        return importlib.import_module(f"spatialcore.{module_name}")
    except ImportError:
        return None


def _module_unavailable(name: str):
    """Create a placeholder that raises helpful error on access."""

    class _UnavailableModule:
        def __init__(self, module_name: str):
            self._name = module_name

        def __getattr__(self, attr: str):
            raise ImportError(
                f"The '{self._name}' module is not available in this version of "
                f"SpatialCore (v{__version__}). Available modules: {_available_modules}. "
                f"Check https://github.com/mcap91/SpatialCore for the latest release."
            )

        def __repr__(self) -> str:
            return f"<UnavailableModule: {self._name}>"

    return _UnavailableModule(name)


# Core module is required - fail loudly if missing
from spatialcore import core

_available_modules.append("core")

# Optional modules - import if available, placeholder if not
# Each module can be shipped independently as it becomes ready

annotation = _try_import("annotation")
if annotation is not None:
    _available_modules.append("annotation")
else:
    annotation = _module_unavailable("annotation")

diffusion = _try_import("diffusion")
if diffusion is not None:
    _available_modules.append("diffusion")
else:
    diffusion = _module_unavailable("diffusion")

nmf = _try_import("nmf")
if nmf is not None:
    _available_modules.append("nmf")
else:
    nmf = _module_unavailable("nmf")

plotting = _try_import("plotting")
if plotting is not None:
    _available_modules.append("plotting")
else:
    plotting = _module_unavailable("plotting")

r_bridge = _try_import("r_bridge")
if r_bridge is not None:
    _available_modules.append("r_bridge")
else:
    r_bridge = _module_unavailable("r_bridge")

spatial = _try_import("spatial")
if spatial is not None:
    _available_modules.append("spatial")
else:
    spatial = _module_unavailable("spatial")

stats = _try_import("stats")
if stats is not None:
    _available_modules.append("stats")
else:
    stats = _module_unavailable("stats")

# Internal modules (only available on dev branch, never published)
internal = _try_import("internal")
if internal is not None:
    _available_modules.append("internal")
else:
    internal = None  # No placeholder - truly internal

__all__ = [
    "__version__",
    "core",
    "annotation",
    "diffusion",
    "nmf",
    "plotting",
    "r_bridge",
    "spatial",
    "stats",
]


def available_modules() -> list[str]:
    """Return list of modules available in this installation."""
    return _available_modules.copy()


def print_info() -> None:
    """Print version and available modules."""
    print(f"SpatialCore v{__version__}")
    print(f"Available modules: {', '.join(_available_modules)}")
