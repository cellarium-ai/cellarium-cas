import logging

NUM_ATTEMPTS_PER_CHUNK_DEFAULT = 7
MAX_NUM_REQUESTS_AT_A_TIME = 8
START_RETRY_DELAY = 5
MAX_RETRY_DELAY = 32
AIOHTTP_TOTAL_TIMEOUT_SECONDS = 750
AIOHTTP_READ_TIMEOUT_SECONDS = 730
MAX_CHUNK_SIZE_SEARCH_METHOD = 500
CELLARIUM_CLOUD_BACKEND_URL = "https://cellarium-cloud-api.cellarium.ai"
LOGGING_LEVEL = logging.INFO
LOGGING_FORMAT = "* [%(asctime)s.%(msecs)03d] %(message)s"
LOGGING_DATE_FORMAT = "%H:%M:%S"


def is_interactive_environment() -> bool:
    """
    Check if the current environment is interactive (e.g. Jupyter notebook, IPython terminal)

    :return: True if the current environment is interactive, False otherwise
    """
    import sys

    import __main__

    # Check if running in an interactive shell (e.g., IPython, Jupyter)
    if hasattr(sys, "ps1"):
        return True
    # Check if running in a context without a '__file__' attribute, typical of interactive environments
    if not hasattr(__main__, "__file__"):
        return True
    # Check specific interactive environment conditions
    try:
        # Check if running in IPython or Jupyter by checking the module name
        from IPython import get_ipython

        shell = get_ipython().__class__.__name__
        if shell in ["ZMQInteractiveShell", "TerminalInteractiveShell"]:
            return True
    except (ModuleNotFoundError, NameError):
        pass
    return False
