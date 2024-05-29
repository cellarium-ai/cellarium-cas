from dataclasses import dataclass


class CASBaseError(Exception):
    """Base class for all CAS exceptions"""


class CASClientError(CASBaseError):
    pass


class HTTPError(CASBaseError):
    pass


class HTTPError5XX(HTTPError):
    pass


class HTTPError404(HTTPError):
    pass


class HTTPError403(HTTPError):
    pass


class HTTPError401(HTTPError):
    pass


class HTTPClientError(HTTPError):
    pass


@dataclass
class DataValidationError(CASBaseError):
    missing_features: int
    extra_features: int


class QuotaExceededError(CASBaseError):
    pass
