from dataclasses import dataclass
from typing import Optional


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
    incompatible_x_type: Optional[str]
    incompatible_total_mrna_umis_type: Optional[str]


class QuotaExceededError(CASBaseError):
    pass


class ClientTooOldError(CASBaseError):
    pass
