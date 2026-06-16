"""
Unit tests for compute_hierarchical_f_measure_from_cm.

Hierarchy used throughout (same as _fixtures.py):
    CL:0000000 (root -- excluded from LABELS and ancestor sets)
        └── CL:0000001 (parent)
                ├── CL:0000002 (child_a)
                └── CL:0000003 (child_b)
                        └── CL:0000004 (grandchild)
"""

import hiclass.metrics as hc_metrics
import numpy as np
import scipy.sparse

from cellarium.cas.benchmarking.confusion_matrix import build_confusion_matrix
from cellarium.cas.benchmarking.hierarchical_f_measure import compute_hierarchical_f_measure_from_cm
from tests.unit.benchmarking._fixtures import LABEL_PATHS, LABELS, ONTOLOGY_CACHE

# ---------------------------------------------------------------------------
# TP / FP / FN (hand-computed)
#
# Example: true=CL:0000002 x2 (correct), true=CL:0000003 x1 (predicted as CL:0000002)
#
#   Correct pair (count=2):
#     T={CL:0000001, CL:0000002}  P={CL:0000001, CL:0000002}
#     T∩P={CL:0000001, CL:0000002} -> hTP=2, hFP=0, hFN=0
#     contribution: hTP+=4, hFP+=0, hFN+=0
#
#   Wrong pair (count=1):
#     T={CL:0000001, CL:0000003}  P={CL:0000001, CL:0000002}
#     T∩P={CL:0000001} -> hTP=1, P-T={CL:0000002} -> hFP=1, T-P={CL:0000003} -> hFN=1
#     contribution: hTP+=1, hFP+=1, hFN+=1
#
#   Totals: hTP=5, hFP=1, hFN=1
# ---------------------------------------------------------------------------


def _make_base_cm():
    """2 correct CL:0000002, 1 wrong CL:0000003 -> CL:0000002."""
    return build_confusion_matrix(
        ["CL:0000002", "CL:0000002", "CL:0000003"],
        ["CL:0000002", "CL:0000002", "CL:0000002"],
        LABELS,
    )


def test_h_tp():
    m = compute_hierarchical_f_measure_from_cm(_make_base_cm(), LABELS, ONTOLOGY_CACHE)
    assert m["h_tp"] == 5


def test_h_fp():
    m = compute_hierarchical_f_measure_from_cm(_make_base_cm(), LABELS, ONTOLOGY_CACHE)
    assert m["h_fp"] == 1


def test_h_fn():
    m = compute_hierarchical_f_measure_from_cm(_make_base_cm(), LABELS, ONTOLOGY_CACHE)
    assert m["h_fn"] == 1


def test_h_tp_fp_fn_perfect():
    # All correct -> hFP = hFN = 0
    cm = build_confusion_matrix(
        ["CL:0000002", "CL:0000002"],
        ["CL:0000002", "CL:0000002"],
        LABELS,
    )
    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
    assert m["h_fp"] == 0
    assert m["h_fn"] == 0


# ---------------------------------------------------------------------------
# Micro metrics
#
# For the base scenario: hTP=5, hFP=1, hFN=1
#   hP = 5/(5+1) = 5/6,  hR = 5/(5+1) = 5/6,  hF1 = 5/6
# ---------------------------------------------------------------------------


def test_h_micro_precision():
    m = compute_hierarchical_f_measure_from_cm(_make_base_cm(), LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_precision_micro"], 5 / 6, rtol=1e-9)


def test_h_micro_recall():
    m = compute_hierarchical_f_measure_from_cm(_make_base_cm(), LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_recall_micro"], 5 / 6, rtol=1e-9)


def test_h_micro_f1():
    m = compute_hierarchical_f_measure_from_cm(_make_base_cm(), LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_f1_micro"], 5 / 6, rtol=1e-9)


def test_h_micro_f1_perfect():
    cm = build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002", "CL:0000003"], LABELS)
    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_f1_micro"], 1.0)


def test_h_micro_vs_hiclass():
    """Validate hierarchical micro metrics against HiClass."""
    y_true = ["CL:0000002", "CL:0000002", "CL:0000003"]
    y_pred = ["CL:0000002", "CL:0000002", "CL:0000002"]
    cm = build_confusion_matrix(y_true, y_pred, LABELS)

    y_true_paths = [LABEL_PATHS[lbl] for lbl in y_true]
    y_pred_paths = [LABEL_PATHS[lbl] for lbl in y_pred]

    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_precision_micro"], hc_metrics.precision(y_true_paths, y_pred_paths), rtol=1e-6)
    np.testing.assert_allclose(m["h_recall_micro"], hc_metrics.recall(y_true_paths, y_pred_paths), rtol=1e-6)
    np.testing.assert_allclose(m["h_f1_micro"], hc_metrics.f1(y_true_paths, y_pred_paths), rtol=1e-6)


def test_h_micro_vs_hiclass_deeper_hierarchy():
    """Validate with grandchild (CL:0000004) predictions using a longer path."""
    y_true = ["CL:0000004", "CL:0000004", "CL:0000002"]
    y_pred = ["CL:0000004", "CL:0000003", "CL:0000001"]
    cm = build_confusion_matrix(y_true, y_pred, LABELS)

    y_true_paths = [LABEL_PATHS[lbl] for lbl in y_true]
    y_pred_paths = [LABEL_PATHS[lbl] for lbl in y_pred]

    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_precision_micro"], hc_metrics.precision(y_true_paths, y_pred_paths), rtol=1e-6)
    np.testing.assert_allclose(m["h_recall_micro"], hc_metrics.recall(y_true_paths, y_pred_paths), rtol=1e-6)
    np.testing.assert_allclose(m["h_f1_micro"], hc_metrics.f1(y_true_paths, y_pred_paths), rtol=1e-6)


# ---------------------------------------------------------------------------
# Macro F1 (per-node OvR, averaged over nodes with nonzero true support)
#
# Base scenario: CL:0000002 x2 correct, CL:0000003 x1 -> CL:0000002
#   Node CL:0000001: tp=3, fp=0, fn=0 -> P=1.0, R=1.0, F1=1.0
#   Node CL:0000002: tp=2, fp=1, fn=0 -> P=2/3, R=1.0, F1=4/5=0.8
#   Node CL:0000003: tp=0, fp=0, fn=1 -> P=0.0, R=0.0, F1=0.0
#   macro F1 = (1.0 + 0.8 + 0.0) / 3 = 0.6
# ---------------------------------------------------------------------------


def test_h_macro_f1_hand_computed():
    m = compute_hierarchical_f_measure_from_cm(_make_base_cm(), LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_f1_macro"], 3 / 5, rtol=1e-9)


def test_h_macro_f1_perfect():
    cm = build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002", "CL:0000003"], LABELS)
    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_f1_macro"], 1.0)


def test_h_macro_uses_only_true_classes():
    # Only CL:0000002 as true label -> macro over one class -> 1.0
    cm = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_f1_macro"], 1.0)


def test_h_macro_symmetric_siblings():
    """
    Each sibling misclassified as the other.
        Cell 1: true=CL:0000002, pred=CL:0000003
        Cell 2: true=CL:0000003, pred=CL:0000002

    Per-node OvR:
        CL:0000001: tp=2, fp=0, fn=0 -> F1=1.0
        CL:0000002: tp=0, fp=1, fn=1 -> F1=0.0
        CL:0000003: tp=0, fp=1, fn=1 -> F1=0.0
    Macro = (1.0 + 0.0 + 0.0) / 3 = 1/3
    """
    cm = build_confusion_matrix(
        ["CL:0000002", "CL:0000003"],
        ["CL:0000003", "CL:0000002"],
        LABELS,
    )
    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_f1_macro"], 1 / 3, rtol=1e-9)


def test_h_macro_f1_zero_division():
    cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE, zero_division=0.0)
    np.testing.assert_allclose(m["h_f1_macro"], 0.0)


# ---------------------------------------------------------------------------
# Weighted metrics
# ---------------------------------------------------------------------------


def test_h_weighted_keys_present():
    m = compute_hierarchical_f_measure_from_cm(_make_base_cm(), LABELS, ONTOLOGY_CACHE)
    for key in ("h_precision_weighted", "h_recall_weighted", "h_f1_weighted"):
        assert key in m, f"Key {key!r} missing from result"
        assert isinstance(m[key], float), f"{key} should be float"
        assert 0.0 <= m[key] <= 1.0, f"{key}={m[key]} out of [0, 1]"


def test_h_weighted_perfect():
    cm = build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002", "CL:0000003"], LABELS)
    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_precision_weighted"], 1.0)
    np.testing.assert_allclose(m["h_recall_weighted"], 1.0)
    np.testing.assert_allclose(m["h_f1_weighted"], 1.0)


def test_h_weighted_empty_cm():
    cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE, zero_division=0.0)
    np.testing.assert_allclose(m["h_precision_weighted"], 0.0)
    np.testing.assert_allclose(m["h_recall_weighted"], 0.0)
    np.testing.assert_allclose(m["h_f1_weighted"], 0.0)


# ---------------------------------------------------------------------------
# Macro precision / recall (hand-computed)
#
# Base scenario:
#   Node CL:0000001: tp=3, fp=0, fn=0 -> P=1.0,  R=1.0
#   Node CL:0000002: tp=2, fp=1, fn=0 -> P=2/3,  R=1.0
#   Node CL:0000003: tp=0, fp=0, fn=1 -> P=0.0,  R=0.0
#   precision_macro = (1.0 + 2/3 + 0.0) / 3 = 5/9
#   recall_macro    = (1.0 + 1.0 + 0.0) / 3 = 2/3
# ---------------------------------------------------------------------------


def test_h_macro_precision_recall_hand_computed():
    m = compute_hierarchical_f_measure_from_cm(_make_base_cm(), LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_precision_macro"], 5 / 9, rtol=1e-9)
    np.testing.assert_allclose(m["h_recall_macro"], 2 / 3, rtol=1e-9)


def test_h_macro_precision_recall_perfect():
    cm = build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002", "CL:0000003"], LABELS)
    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
    np.testing.assert_allclose(m["h_precision_macro"], 1.0)
    np.testing.assert_allclose(m["h_recall_macro"], 1.0)


def test_h_macro_precision_recall_empty_cm():
    cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
    m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE, zero_division=0.0)
    np.testing.assert_allclose(m["h_precision_macro"], 0.0)
    np.testing.assert_allclose(m["h_recall_macro"], 0.0)
