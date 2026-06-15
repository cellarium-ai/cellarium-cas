"""
Hierarchical F-measure metrics computed from a confusion matrix.

Implements the Kiritchenko et al. hierarchical evaluation approach: each cell
in the confusion matrix contributes hierarchical TP, FP, and FN counts based on
the overlap between the ontology ancestor sets of the true and predicted labels.

    For each nonzero cell M[true, pred] with count *c*:

        T = ancestors(true_label)   # self-inclusive, root-exclusive
        P = ancestors(pred_label)

        hTP += c * |T ∩ P|
        hFP += c * |P − T|
        hFN += c * |T − P|

Micro metrics aggregate these counts globally.  Macro metrics compute per-node
(OvR) F1 scores across all ontology ancestor nodes with nonzero true support,
then average them.
"""

import typing as t
from collections import defaultdict

import numpy as np
import scipy.sparse

from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache


def compute_hierarchical_f_measure_from_cm(
    cm: scipy.sparse.csr_matrix,
    label_order: t.List[str],
    ontology_cache: CellOntologyCache,
    zero_division: float = 0.0,
) -> t.Dict[str, t.Any]:
    """
    Compute hierarchical F-measure metrics from a confusion matrix.

    :param cm: Square sparse confusion matrix of shape ``(N, N)`` where N =
        ``len(label_order)``.  Row *i* = true label at ``label_order[i]``;
        column *j* = predicted label.
    :param label_order: Ordered list of labels corresponding to the CM rows/cols.
        Must match the ``label_order`` used when building the CM.
    :param ontology_cache: ``CellOntologyCache`` instance used to look up ancestor sets
        (self-inclusive, root-exclusive) for each label.
    :param zero_division: Value to return when a denominator is zero (default ``0.0``).
    :returns: Dict with keys:

        - ``h_tp`` — global hierarchical true positives (``int``).
        - ``h_fp`` — global hierarchical false positives (``int``).
        - ``h_fn`` — global hierarchical false negatives (``int``).
        - ``h_precision_micro`` — micro hierarchical precision.
        - ``h_recall_micro`` — micro hierarchical recall.
        - ``h_f1_micro`` — micro hierarchical F1.
        - ``h_f1_macro`` — macro hierarchical F1 over ontology nodes with nonzero true support (OvR).
    """
    cm_coo = cm.tocoo()

    h_tp = 0.0
    h_fp = 0.0
    h_fn = 0.0

    node_h_tp: t.Dict[str, float] = defaultdict(float)
    node_h_fp: t.Dict[str, float] = defaultdict(float)
    node_h_fn: t.Dict[str, float] = defaultdict(float)

    for idx in range(len(cm_coo.data)):
        count = float(cm_coo.data[idx])
        if count == 0:
            continue
        i = int(cm_coo.row[idx])
        j = int(cm_coo.col[idx])
        true_label = label_order[i]
        pred_label = label_order[j]

        true_ancestors_set = ontology_cache.get_ancestors(true_label, remove_root=True)
        pred_ancestors_set = ontology_cache.get_ancestors(pred_label, remove_root=True)

        tp_set = true_ancestors_set & pred_ancestors_set
        fp_set = pred_ancestors_set - true_ancestors_set
        fn_set = true_ancestors_set - pred_ancestors_set

        h_tp += count * len(tp_set)
        h_fp += count * len(fp_set)
        h_fn += count * len(fn_set)

        for node in tp_set:
            node_h_tp[node] += count
        for node in fp_set:
            node_h_fp[node] += count
        for node in fn_set:
            node_h_fn[node] += count

    # --- Micro metrics ---
    h_precision_micro = _safe_div(h_tp, h_tp + h_fp, zero_division)
    h_recall_micro = _safe_div(h_tp, h_tp + h_fn, zero_division)
    h_f1_micro = _safe_f1(h_precision_micro, h_recall_micro, zero_division)

    # --- Macro metrics: OvR per ontology node, averaged over nodes with nonzero true support ---
    node_f1s: t.List[float] = []
    for node in set(node_h_tp.keys()) | set(node_h_fn.keys()):
        tp = node_h_tp[node]
        fn = node_h_fn[node]
        if tp + fn == 0:
            continue
        fp = node_h_fp[node]
        p = _safe_div(tp, tp + fp, zero_division)
        r = _safe_div(tp, tp + fn, zero_division)
        node_f1s.append(_safe_f1(p, r, zero_division))

    h_f1_macro = float(np.mean(node_f1s)) if node_f1s else zero_division

    return {
        "h_tp": int(h_tp),
        "h_fp": int(h_fp),
        "h_fn": int(h_fn),
        "h_precision_micro": float(h_precision_micro),
        "h_recall_micro": float(h_recall_micro),
        "h_f1_micro": float(h_f1_micro),
        "h_f1_macro": float(h_f1_macro),
    }


def _safe_div(numerator: float, denominator: float, zero_division: float) -> float:
    return numerator / denominator if denominator > 0 else zero_division


def _safe_f1(precision: float, recall: float, zero_division: float) -> float:
    denom = precision + recall
    return 2.0 * precision * recall / denom if denom > 0 else zero_division
