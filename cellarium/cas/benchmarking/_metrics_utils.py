"""Shared arithmetic helpers for F-measure and hierarchical F-measure modules."""


def _safe_div(numerator: float, denominator: float, zero_division: float) -> float:
    return numerator / denominator if denominator > 0 else zero_division


def _safe_f1(precision: float, recall: float, zero_division: float) -> float:
    denom = precision + recall
    return 2.0 * precision * recall / denom if denom > 0 else zero_division
