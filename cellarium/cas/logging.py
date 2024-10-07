import logging

from cellarium.cas import settings

# Create a custom logger for the package
logger = logging.getLogger("cellarium.cas")
logger.setLevel(settings.LOGGING_LEVEL)

# Create a handler
handler = logging.StreamHandler()

# Create a formatter and set it for the handler
formatter = logging.Formatter(fmt=settings.LOGGING_FORMAT, datefmt=settings.LOGGING_DATE_FORMAT)
handler.setFormatter(formatter)

# Add the handler to the logger
logger.addHandler(handler)
