import os
from sys import version_info as python_version_info


def get_version() -> str:
    """
    Get the version of the application based on the git tag.

    :return: The version of the application.
    """

    if python_version_info.major == 3 and python_version_info.minor >= 8:
        from importlib.metadata import PackageNotFoundError, version

        try:
            return version("cellarium-cas")
        except PackageNotFoundError:
            return os.environ.get("CAS_VERSION", "0.0.1")
    elif python_version_info.major == 3 and python_version_info.minor == 7:
        from importlib_metadata import PackageNotFoundError, version

        try:
            return version("cellarium-cas")
        except PackageNotFoundError:
            return os.environ.get("CAS_VERSION", "0.0.1")
    else:
        raise Exception("Unsupported Python version. Please use Python 3.7 or later.")
