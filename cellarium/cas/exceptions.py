from dataclasses import dataclass


class CASBaseError(Exception):
    """Base class for all CAS exceptions"""


class CASClientError(CASBaseError):
    pass


class HTTPError(CASBaseError):
    pass


class HTTPError500(HTTPError):
    pass


class HTTPError403(HTTPError):
    pass


class HTTPError401(HTTPError):
    pass


@dataclass
class DataValidationError(CASBaseError):
    missing_features: int
    extra_features: int
