NUM_ATTEMPTS_PER_CHUNK_DEFAULT = 7
MAX_NUM_REQUESTS_AT_A_TIME = 25
START_RETRY_DELAY = 5
MAX_RETRY_DELAY = 32
AIOHTTP_TOTAL_TIMEOUT_SECONDS = 750
AIOHTTP_READ_TIMEOUT_SECONDS = 730
MAX_CHUNK_SIZE_SEARCH_METHOD = 500


def is_interactive_environment() -> bool:
    """
    Check if the current environment is interactive (e.g. Jupyter notebook, IPython terminal)

    :return: True if the current environment is interactive, False otherwise
    """
    try:
        from IPython import get_ipython

        if "IPKernelApp" not in get_ipython().config:
            # Running in a non-interactive environment.
            return False
        else:
            if "terminal" in get_ipython().config["IPKernelApp"]["connection_file"]:
                # Running in an IPython terminal
                return True
            else:
                # Running in a Jupyter environment.
                return True
    except (ImportError, AttributeError):
        # Running in a non-interactive environment.
        return False
