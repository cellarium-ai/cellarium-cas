"""
Unit tests for compute_f_measure_from_cm.
"""

import numpy as np
import scipy.sparse
from sklearn.metrics import f1_score, precision_score, recall_score

from cellarium.cas.benchmarking.confusion_matrix import build_confusion_matrix
from cellarium.cas.benchmarking.f_measure import compute_f_measure_from_cm
from tests.unit.benchmarking._fixtures import LABELS

# ---------------------------------------------------------------------------
# compute_f_measure_from_cm
# ---------------------------------------------------------------------------


def test_f_measure_perfect_predictions():
    y_true = ["CL:0000002", "CL:0000002", "CL:0000003"]
    y_pred = ["CL:0000002", "CL:0000002", "CL:0000003"]
    cm = build_confusion_matrix(y_true, y_pred, LABELS)
    m = compute_f_measure_from_cm(cm)
    assert m["tp"] == 3
    assert m["fp"] == 0
    assert m["fn"] == 0
    np.testing.assert_allclose(m["precision_micro"], 1.0)
    np.testing.assert_allclose(m["recall_micro"], 1.0)
    np.testing.assert_allclose(m["f1_micro"], 1.0)
    np.testing.assert_allclose(m["f1_macro"], 1.0)


def test_f_measure_all_wrong():
    y_true = ["CL:0000002", "CL:0000002"]
    y_pred = ["CL:0000003", "CL:0000003"]
    cm = build_confusion_matrix(y_true, y_pred, LABELS)
    m = compute_f_measure_from_cm(cm)
    assert m["tp"] == 0
    assert m["fp"] == 2
    assert m["fn"] == 2
    np.testing.assert_allclose(m["f1_micro"], 0.0)


def test_f_measure_partial_predictions():
    # 2 correct, 1 wrong
    y_true = ["CL:0000002", "CL:0000002", "CL:0000003"]
    y_pred = ["CL:0000002", "CL:0000002", "CL:0000002"]
    cm = build_confusion_matrix(y_true, y_pred, LABELS)
    m = compute_f_measure_from_cm(cm)
    assert m["tp"] == 2
    assert m["fp"] == 1
    assert m["fn"] == 1
    np.testing.assert_allclose(m["precision_micro"], 2 / 3)
    np.testing.assert_allclose(m["recall_micro"], 2 / 3)
    np.testing.assert_allclose(m["f1_micro"], 2 / 3)


def test_f_measure_macro_excludes_zero_support_classes():
    # Only CL:0000002 appears; macro should be over one class only
    cm = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
    m = compute_f_measure_from_cm(cm)
    np.testing.assert_allclose(m["f1_macro"], 1.0)


def test_f_measure_zero_division_default_is_zero():
    cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
    m = compute_f_measure_from_cm(cm)
    np.testing.assert_allclose(m["f1_micro"], 0.0)
    np.testing.assert_allclose(m["f1_macro"], 0.0)


def test_f_measure_zero_division_custom_value():
    cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
    m = compute_f_measure_from_cm(cm, zero_division=1.0)
    np.testing.assert_allclose(m["f1_micro"], 1.0)


def test_f_measure_weighted_keys_present_and_perfect():
    y_true = ["CL:0000002", "CL:0000002", "CL:0000003"]
    y_pred = ["CL:0000002", "CL:0000002", "CL:0000003"]
    cm = build_confusion_matrix(y_true, y_pred, LABELS)
    m = compute_f_measure_from_cm(cm)
    np.testing.assert_allclose(m["precision_weighted"], 1.0)
    np.testing.assert_allclose(m["recall_weighted"], 1.0)
    np.testing.assert_allclose(m["f1_weighted"], 1.0)


def test_f_measure_weighted_matches_sklearn():
    y_true = ["CL:0000002", "CL:0000002", "CL:0000003", "CL:0000004"]
    y_pred = ["CL:0000002", "CL:0000003", "CL:0000003", "CL:0000004"]
    cm = build_confusion_matrix(y_true, y_pred, LABELS)
    m = compute_f_measure_from_cm(cm)
    np.testing.assert_allclose(
        m["precision_weighted"],
        precision_score(y_true, y_pred, labels=LABELS, average="weighted", zero_division=0),
        rtol=1e-6,
    )
    np.testing.assert_allclose(
        m["recall_weighted"],
        recall_score(y_true, y_pred, labels=LABELS, average="weighted", zero_division=0),
        rtol=1e-6,
    )
    np.testing.assert_allclose(
        m["f1_weighted"],
        f1_score(y_true, y_pred, labels=LABELS, average="weighted", zero_division=0),
        rtol=1e-6,
    )


def test_f_measure_weighted_zero_division_default():
    cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
    m = compute_f_measure_from_cm(cm)
    np.testing.assert_allclose(m["precision_weighted"], 0.0)
    np.testing.assert_allclose(m["recall_weighted"], 0.0)
    np.testing.assert_allclose(m["f1_weighted"], 0.0)


def test_f_measure_weighted_zero_division_custom():
    cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
    m = compute_f_measure_from_cm(cm, zero_division=1.0)
    np.testing.assert_allclose(m["precision_weighted"], 1.0)
    np.testing.assert_allclose(m["recall_weighted"], 1.0)
    np.testing.assert_allclose(m["f1_weighted"], 1.0)


def test_f_measure_macro_precision_recall_perfect():
    y_true = ["CL:0000002", "CL:0000002", "CL:0000003"]
    y_pred = ["CL:0000002", "CL:0000002", "CL:0000003"]
    cm = build_confusion_matrix(y_true, y_pred, LABELS)
    m = compute_f_measure_from_cm(cm)
    np.testing.assert_allclose(m["precision_macro"], 1.0)
    np.testing.assert_allclose(m["recall_macro"], 1.0)


def test_f_measure_macro_precision_recall_matches_sklearn():
    # Our macro averages over nonzero-support classes only; pass those same labels to sklearn.
    y_true = ["CL:0000002", "CL:0000002", "CL:0000003", "CL:0000004"]
    y_pred = ["CL:0000002", "CL:0000003", "CL:0000003", "CL:0000004"]
    active_labels = sorted(set(y_true) | set(y_pred))
    cm = build_confusion_matrix(y_true, y_pred, LABELS)
    m = compute_f_measure_from_cm(cm)
    np.testing.assert_allclose(
        m["precision_macro"],
        precision_score(y_true, y_pred, labels=active_labels, average="macro", zero_division=0),
        rtol=1e-6,
    )
    np.testing.assert_allclose(
        m["recall_macro"],
        recall_score(y_true, y_pred, labels=active_labels, average="macro", zero_division=0),
        rtol=1e-6,
    )


def test_f_measure_macro_precision_recall_zero_division_default():
    cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
    m = compute_f_measure_from_cm(cm)
    np.testing.assert_allclose(m["precision_macro"], 0.0)
    np.testing.assert_allclose(m["recall_macro"], 0.0)


def test_f_measure_macro_precision_recall_zero_division_custom():
    cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
    m = compute_f_measure_from_cm(cm, zero_division=1.0)
    np.testing.assert_allclose(m["precision_macro"], 1.0)
    np.testing.assert_allclose(m["recall_macro"], 1.0)
