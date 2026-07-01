from cellarium.cas.client import CASClient

from . import constants, exceptions, logging, postprocessing, preprocessing, service, settings, version, visualization

try:
    from . import benchmarking
except ImportError:
    benchmarking = None  # type: ignore[assignment]

__version__ = version.get_version()

__all__ = [
    "CASClient",
    "benchmarking",
    "preprocessing",
    "postprocessing",
    "visualization",
    "settings",
    "constants",
    "exceptions",
    "service",
    "settings",
    "version",
    "logging",
]
