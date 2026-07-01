"""
Standard (non-hierarchical) F-measure metrics computed from a confusion matrix.

Micro precision, recall, and F1 treat the problem globally: TP is the trace of the
confusion matrix, FP = FN = total − trace.  For single-label multiclass classification
this means micro precision = recall = F1 = accuracy.

Macro metrics are computed over labels present in either y_true or y_pred.
Weighted metrics use true-label support as weights.
"""

import typing as t

import numpy as np
import scipy.sparse

from ._metrics_utils import _safe_div, _safe_f1


def compute_f_measure_from_cm(
    cm: scipy.sparse.csr_matrix,
    zero_division: float = 0.0,
) -> t.Dict[str, t.Any]:
    """
    Compute standard F-measure metrics from a (possibly sparse) confusion matrix.

    :param cm: Square sparse confusion matrix, rows = true labels, cols = predicted labels.
    :param zero_division: Value to return when a denominator is zero (default ``0.0``).
        Mirrors the ``zero_division`` parameter in :mod:`sklearn.metrics`.
    :returns: Dict with keys:

        - ``tp`` — global true positives (``int``).
        - ``fp`` — global false positives (``int``).
        - ``fn`` — global false negatives (``int``).
        - ``precision_micro`` — micro precision.
        - ``recall_micro`` — micro recall.
        - ``f1_micro`` — micro F1.
        - ``f1_macro`` — macro F1 over nonzero-support classes.
        - ``precision_macro`` — macro precision over nonzero-support classes.
        - ``recall_macro`` — macro recall over nonzero-support classes.
        - ``precision_weighted`` — support-weighted precision over nonzero-support classes.
        - ``recall_weighted`` — support-weighted recall over nonzero-support classes.
        - ``f1_weighted`` — support-weighted F1 over nonzero-support classes.
    """
    row_sums = np.asarray(cm.sum(axis=1)).ravel().astype(np.int64)
    col_sums = np.asarray(cm.sum(axis=0)).ravel().astype(np.int64)
    diag = cm.diagonal().astype(np.int64)

    tp = int(diag.sum())
    fp = int((col_sums - diag).sum())
    fn = int((row_sums - diag).sum())

    precision_micro = _safe_div(tp, tp + fp, zero_division)
    recall_micro = _safe_div(tp, tp + fn, zero_division)
    f1_micro = _safe_f1(precision_micro, recall_micro, zero_division)

    nonzero_mask = (row_sums > 0) | (col_sums > 0)
    class_tp = diag[nonzero_mask].astype(np.float64)
    class_fp = (col_sums[nonzero_mask] - diag[nonzero_mask]).astype(np.float64)
    class_fn = (row_sums[nonzero_mask] - diag[nonzero_mask]).astype(np.float64)

    denom_p = class_tp + class_fp
    denom_r = class_tp + class_fn
    class_p = np.where(denom_p > 0, class_tp / denom_p, zero_division)
    class_r = np.where(denom_r > 0, class_tp / denom_r, zero_division)
    denom_f1 = class_p + class_r
    class_f1 = np.where(denom_f1 > 0, 2.0 * class_p * class_r / denom_f1, zero_division)

    f1_macro = float(class_f1.mean()) if class_f1.size > 0 else zero_division
    precision_macro = float(class_p.mean()) if class_p.size > 0 else zero_division
    recall_macro = float(class_r.mean()) if class_r.size > 0 else zero_division

    support = row_sums[nonzero_mask].astype(np.float64)
    total_support = support.sum()
    if total_support > 0:
        precision_weighted = float((class_p * support).sum() / total_support)
        recall_weighted = float((class_r * support).sum() / total_support)
        f1_weighted = float((class_f1 * support).sum() / total_support)
    else:
        precision_weighted = recall_weighted = f1_weighted = zero_division

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision_micro": float(precision_micro),
        "recall_micro": float(recall_micro),
        "f1_micro": float(f1_micro),
        "f1_macro": float(f1_macro),
        "precision_macro": float(precision_macro),
        "recall_macro": float(recall_macro),
        "precision_weighted": float(precision_weighted),
        "recall_weighted": float(recall_weighted),
        "f1_weighted": float(f1_weighted),
    }
