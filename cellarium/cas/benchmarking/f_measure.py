"""
Standard (non-hierarchical) F-measure metrics computed from a confusion matrix.

Micro precision, recall, and F1 treat the problem globally: TP is the trace of the
confusion matrix, FP = FN = total − trace.  For single-label multiclass classification
this means micro precision = recall = F1 = accuracy.

Macro F1 is the unweighted mean of per-class F1 scores, computed only over classes
with nonzero support (i.e. appearing as a true label **or** as a predicted label in
the matrix).
"""

import typing as t

import numpy as np
import scipy.sparse


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
    """
    dense = cm.toarray().astype(np.int64)
    total = int(dense.sum())
    tp = int(np.trace(dense))
    fp = total - tp
    fn = total - tp

    precision_micro = _safe_div(tp, tp + fp, zero_division)
    recall_micro = _safe_div(tp, tp + fn, zero_division)
    f1_micro = _safe_f1(precision_micro, recall_micro, zero_division)

    row_sums = dense.sum(axis=1)
    col_sums = dense.sum(axis=0)
    diag = np.diag(dense)

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

    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision_micro": float(precision_micro),
        "recall_micro": float(recall_micro),
        "f1_micro": float(f1_micro),
        "f1_macro": float(f1_macro),
    }


def _safe_div(numerator: float, denominator: float, zero_division: float) -> float:
    return numerator / denominator if denominator > 0 else zero_division


def _safe_f1(precision: float, recall: float, zero_division: float) -> float:
    denom = precision + recall
    return 2.0 * precision * recall / denom if denom > 0 else zero_division
