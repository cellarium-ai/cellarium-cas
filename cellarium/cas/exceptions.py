from dataclasses import dataclass


class CASBaseError(Exception):
    """Base class for all CAS exceptions"""

    pass


class CASClientError(CASBaseError):
    pass


class HTTPBaseError(CASBaseError):
    pass


class HTTPError500(HTTPBaseError):
    pass


class HTTPError403(HTTPBaseError):
    pass


class HTTPError401(HTTPBaseError):
    pass


@dataclass
class DataValidationError(CASBaseError):
    missing_features: int
    extra_features: int
