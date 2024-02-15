from cellarium.cas.settings.base import *  # noqa
from cellarium.cas.version import get_version_environment

VERSION_ENVIRONMENT = get_version_environment()

if VERSION_ENVIRONMENT == "development" or VERSION_ENVIRONMENT == "test":
    from cellarium.cas.settings.development import *  # noqa
else:
    from cellarium.cas.settings.production import *  # noqa
