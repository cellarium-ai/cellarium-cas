from dataclasses import dataclass


class HTTPBaseError(Exception):
    pass


class HTTPError500(HTTPBaseError):
    pass


class HTTPError403(HTTPBaseError):
    pass


class HTTPError401(HTTPBaseError):
    pass


@dataclass
class DataValidationError(Exception):
    missing_features: int
    extra_features: int
