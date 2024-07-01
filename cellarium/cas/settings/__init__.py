from ..version import get_version_environment
from .base import *  # noqa

VERSION_ENVIRONMENT = get_version_environment()

if VERSION_ENVIRONMENT == "development" or VERSION_ENVIRONMENT == "test":
    from .development import *  # noqa
else:
    from .production import *  # noqa
