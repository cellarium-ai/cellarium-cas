"""
Unit tests for the confusion-matrix-based F-measure and hierarchical F-measure metrics.

Hierarchy used throughout:

    CL:0000000 (root -- excluded from all ancestor sets)
        └── CL:0000001 (parent)
                ├── CL:0000002 (child_a)
                └── CL:0000003 (child_b)
                        └── CL:0000004 (grandchild)

Ancestor sets (self-inclusive, root-exclusive):
    CL:0000001 -> {CL:0000001}
    CL:0000002 -> {CL:0000001, CL:0000002}
    CL:0000003 -> {CL:0000001, CL:0000003}
    CL:0000004 -> {CL:0000001, CL:0000003, CL:0000004}
"""

import typing as t

import hiclass.metrics as hc_metrics
import numpy as np
import pytest
import scipy.sparse
from sklearn.metrics import f1_score, precision_score, recall_score

from cellarium.cas.benchmarking.confusion_matrix import (
    aggregate_confusion_matrices,
    build_confusion_matrix,
    load_confusion_matrix,
    save_confusion_matrix,
)
from cellarium.cas.benchmarking.f_measure import compute_f_measure_from_cm
from cellarium.cas.benchmarking.hierarchical_f_measure import compute_hierarchical_f_measure_from_cm
from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache

# ---------------------------------------------------------------------------
# Fixtures / shared data
# ---------------------------------------------------------------------------

# Label universe (excluding root so it never appears in any set)
LABELS: t.List[str] = ["CL:0000001", "CL:0000002", "CL:0000003", "CL:0000004"]

# Ancestor sets: self-inclusive, root-exclusive
ANCESTORS: t.Dict[str, t.FrozenSet[str]] = {
    "CL:0000001": frozenset({"CL:0000001"}),
    "CL:0000002": frozenset({"CL:0000001", "CL:0000002"}),
    "CL:0000003": frozenset({"CL:0000001", "CL:0000003"}),
    "CL:0000004": frozenset({"CL:0000001", "CL:0000003", "CL:0000004"}),
}

# Label paths for HiClass (each path = ancestor list from shallowest non-root to leaf)
LABEL_PATHS: t.Dict[str, t.List[str]] = {
    "CL:0000001": ["CL:0000001"],
    "CL:0000002": ["CL:0000001", "CL:0000002"],
    "CL:0000003": ["CL:0000001", "CL:0000003"],
    "CL:0000004": ["CL:0000001", "CL:0000003", "CL:0000004"],
}

MOCK_ONTOLOGY_RESOURCE: t.Dict[str, t.Any] = {
    "cl_names": ["CL:0000000", "CL:0000001", "CL:0000002", "CL:0000003", "CL:0000004"],
    "cell_ontology_term_id_to_cell_type": {
        "CL:0000000": "cell",
        "CL:0000001": "parent",
        "CL:0000002": "child_a",
        "CL:0000003": "child_b",
        "CL:0000004": "grandchild",
    },
    "children_dictionary": {
        "CL:0000000": ["CL:0000001"],
        "CL:0000001": ["CL:0000002", "CL:0000003"],
        "CL:0000003": ["CL:0000004"],
    },
    "shortest_path_lengths_from_cell_root": {
        "CL:0000000": 0,
        "CL:0000001": 1,
        "CL:0000002": 2,
        "CL:0000003": 2,
        "CL:0000004": 3,
    },
    "longest_path_lengths_from_cell_root": {
        "CL:0000000": 0,
        "CL:0000001": 1,
        "CL:0000002": 2,
        "CL:0000003": 2,
        "CL:0000004": 3,
    },
}

ONTOLOGY_CACHE = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)

# ---------------------------------------------------------------------------
# build_confusion_matrix
# ---------------------------------------------------------------------------


class TestBuildConfusionMatrix:
    def test_returns_sparse_matrix(self):
        cm = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
        assert scipy.sparse.issparse(cm)

    def test_shape_equals_label_count(self):
        cm = build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002", "CL:0000002"], LABELS)
        assert cm.shape == (len(LABELS), len(LABELS))

    def test_correct_prediction_on_diagonal(self):
        # CL:0000002 (index 1) predicted correctly
        cm = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
        dense = cm.toarray()
        assert dense[1, 1] == 1  # row=true=CL:0000002, col=pred=CL:0000002

    def test_off_diagonal_for_wrong_prediction(self):
        # true=CL:0000003 (idx 2), pred=CL:0000002 (idx 1)
        cm = build_confusion_matrix(["CL:0000003"], ["CL:0000002"], LABELS)
        dense = cm.toarray()
        assert dense[2, 1] == 1
        assert dense[1, 1] == 0

    def test_multiple_cells_correct_counts(self):
        y_true = ["CL:0000002", "CL:0000002", "CL:0000003"]
        y_pred = ["CL:0000002", "CL:0000002", "CL:0000002"]
        cm = build_confusion_matrix(y_true, y_pred, LABELS)
        dense = cm.toarray()
        assert dense[1, 1] == 2  # 2 correct for CL:0000002
        assert dense[2, 1] == 1  # 1 wrong: CL:0000003 -> CL:0000002

    def test_raises_on_unmapped_gt_label(self):
        with pytest.raises(ValueError, match="not in the ontology"):
            build_confusion_matrix(["CL:9999999"], ["CL:0000002"], LABELS)

    def test_raises_on_unmapped_pred_label(self):
        with pytest.raises(ValueError, match="not in the ontology"):
            build_confusion_matrix(["CL:0000002"], ["CL:9999999"], LABELS)

    def test_raises_on_length_mismatch(self):
        with pytest.raises(ValueError, match="same length"):
            build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002"], LABELS)

    def test_all_labels_in_universe_produce_correct_shape(self):
        # Full-universe label order including root
        full_labels = MOCK_ONTOLOGY_RESOURCE["cl_names"]
        cm = build_confusion_matrix(["CL:0000001"], ["CL:0000001"], full_labels)
        assert cm.shape == (len(full_labels), len(full_labels))


# ---------------------------------------------------------------------------
# aggregate_confusion_matrices
# ---------------------------------------------------------------------------


class TestAggregateConfusionMatrices:
    def test_sum_of_two_matrices(self):
        cm1 = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
        cm2 = build_confusion_matrix(["CL:0000003"], ["CL:0000002"], LABELS)
        agg = aggregate_confusion_matrices([cm1, cm2])
        dense = agg.toarray()
        assert dense[1, 1] == 1  # one correct CL:0000002
        assert dense[2, 1] == 1  # one wrong: CL:0000003 -> CL:0000002

    def test_sum_preserves_shape(self):
        cm1 = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
        cm2 = build_confusion_matrix(["CL:0000001"], ["CL:0000001"], LABELS)
        agg = aggregate_confusion_matrices([cm1, cm2])
        assert agg.shape == cm1.shape

    def test_single_matrix_unchanged(self):
        cm = build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002", "CL:0000003"], LABELS)
        agg = aggregate_confusion_matrices([cm])
        np.testing.assert_array_equal(agg.toarray(), cm.toarray())

    def test_raises_on_empty_list(self):
        with pytest.raises(ValueError, match="empty"):
            aggregate_confusion_matrices([])

    def test_raises_on_shape_mismatch(self):
        cm1 = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
        cm2 = scipy.sparse.csr_matrix(np.zeros((3, 3), dtype=np.int64))
        with pytest.raises(ValueError, match="shape"):
            aggregate_confusion_matrices([cm1, cm2])


# ---------------------------------------------------------------------------
# save / load confusion matrix (round-trip)
# ---------------------------------------------------------------------------


class TestSaveLoadConfusionMatrix:
    def test_roundtrip(self, tmp_path):
        cm = build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002", "CL:0000002"], LABELS)
        meta = {
            "model_name": "test_model",
            "test_sample": "/data/sample.h5ad",
            "label_order": LABELS,
            "matrix_shape": list(cm.shape),
        }
        save_confusion_matrix(cm, meta, tmp_path)
        loaded_cm, loaded_meta = load_confusion_matrix(tmp_path)

        np.testing.assert_array_equal(loaded_cm.toarray(), cm.toarray())
        assert loaded_meta["model_name"] == "test_model"
        assert loaded_meta["label_order"] == LABELS

    def test_load_missing_matrix_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError, match="confusion_matrix.npz"):
            load_confusion_matrix(tmp_path)


# ---------------------------------------------------------------------------
# compute_f_measure_from_cm
# ---------------------------------------------------------------------------


class TestFMeasureFromCM:
    def test_perfect_predictions(self):
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

    def test_all_wrong(self):
        y_true = ["CL:0000002", "CL:0000002"]
        y_pred = ["CL:0000003", "CL:0000003"]
        cm = build_confusion_matrix(y_true, y_pred, LABELS)
        m = compute_f_measure_from_cm(cm)
        assert m["tp"] == 0
        assert m["fp"] == 2
        assert m["fn"] == 2
        np.testing.assert_allclose(m["f1_micro"], 0.0)

    def test_partial_predictions(self):
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

    def test_macro_excludes_zero_support_classes(self):
        # Only CL:0000002 appears; macro should be over one class only
        cm = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
        m = compute_f_measure_from_cm(cm)
        np.testing.assert_allclose(m["f1_macro"], 1.0)

    def test_zero_division_default_is_zero(self):
        # All-zero CM (no predictions)
        cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
        m = compute_f_measure_from_cm(cm)
        np.testing.assert_allclose(m["f1_micro"], 0.0)
        np.testing.assert_allclose(m["f1_macro"], 0.0)

    def test_zero_division_custom_value(self):
        cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
        m = compute_f_measure_from_cm(cm, zero_division=1.0)
        np.testing.assert_allclose(m["f1_micro"], 1.0)

    def test_weighted_keys_present_and_perfect(self):
        y_true = ["CL:0000002", "CL:0000002", "CL:0000003"]
        y_pred = ["CL:0000002", "CL:0000002", "CL:0000003"]
        cm = build_confusion_matrix(y_true, y_pred, LABELS)
        m = compute_f_measure_from_cm(cm)
        np.testing.assert_allclose(m["precision_weighted"], 1.0)
        np.testing.assert_allclose(m["recall_weighted"], 1.0)
        np.testing.assert_allclose(m["f1_weighted"], 1.0)

    def test_weighted_matches_sklearn(self):
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

    def test_weighted_zero_division(self):
        cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
        m = compute_f_measure_from_cm(cm)
        np.testing.assert_allclose(m["precision_weighted"], 0.0)
        np.testing.assert_allclose(m["recall_weighted"], 0.0)
        np.testing.assert_allclose(m["f1_weighted"], 0.0)
        m1 = compute_f_measure_from_cm(cm, zero_division=1.0)
        np.testing.assert_allclose(m1["precision_weighted"], 1.0)
        np.testing.assert_allclose(m1["recall_weighted"], 1.0)
        np.testing.assert_allclose(m1["f1_weighted"], 1.0)

    def test_macro_precision_recall_perfect(self):
        y_true = ["CL:0000002", "CL:0000002", "CL:0000003"]
        y_pred = ["CL:0000002", "CL:0000002", "CL:0000003"]
        cm = build_confusion_matrix(y_true, y_pred, LABELS)
        m = compute_f_measure_from_cm(cm)
        np.testing.assert_allclose(m["precision_macro"], 1.0)
        np.testing.assert_allclose(m["recall_macro"], 1.0)

    def test_macro_precision_recall_matches_sklearn(self):
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

    def test_macro_precision_recall_zero_division(self):
        cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
        m = compute_f_measure_from_cm(cm)
        np.testing.assert_allclose(m["precision_macro"], 0.0)
        np.testing.assert_allclose(m["recall_macro"], 0.0)
        m1 = compute_f_measure_from_cm(cm, zero_division=1.0)
        np.testing.assert_allclose(m1["precision_macro"], 1.0)
        np.testing.assert_allclose(m1["recall_macro"], 1.0)


# ---------------------------------------------------------------------------
# compute_hierarchical_f_measure_from_cm — TP/FP/FN
# ---------------------------------------------------------------------------


class TestHierarchicalTPFPFN:
    """
    Hand-computed ground truth for the hierarchical TP/FP/FN counts.

    Example used:
        Cells: true=CL:0000002 x2 (predicted correctly), true=CL:0000003 x1 (predicted as CL:0000002)

        Correct pair (count=2):
            T = {CL:0000001, CL:0000002}   P = {CL:0000001, CL:0000002}
            T∩P = {CL:0000001, CL:0000002} -> hTP_pair = 2
            P-T = {}                        -> hFP_pair = 0
            T-P = {}                        -> hFN_pair = 0
            contribution: hTP += 4, hFP += 0, hFN += 0

        Wrong pair (count=1):
            T = {CL:0000001, CL:0000003}   P = {CL:0000001, CL:0000002}
            T∩P = {CL:0000001}             -> hTP_pair = 1
            P-T = {CL:0000002}             -> hFP_pair = 1
            T-P = {CL:0000003}             -> hFN_pair = 1
            contribution: hTP += 1, hFP += 1, hFN += 1

        Totals: hTP=5, hFP=1, hFN=1
    """

    def _make_cm(self):
        return build_confusion_matrix(
            ["CL:0000002", "CL:0000002", "CL:0000003"],
            ["CL:0000002", "CL:0000002", "CL:0000002"],
            LABELS,
        )

    def test_h_tp(self):
        m = compute_hierarchical_f_measure_from_cm(self._make_cm(), LABELS, ONTOLOGY_CACHE)
        assert m["h_tp"] == 5

    def test_h_fp(self):
        m = compute_hierarchical_f_measure_from_cm(self._make_cm(), LABELS, ONTOLOGY_CACHE)
        assert m["h_fp"] == 1

    def test_h_fn(self):
        m = compute_hierarchical_f_measure_from_cm(self._make_cm(), LABELS, ONTOLOGY_CACHE)
        assert m["h_fn"] == 1

    def test_perfect_correct_case(self):
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
# compute_hierarchical_f_measure_from_cm — micro metrics
# ---------------------------------------------------------------------------


class TestHierarchicalMicro:
    """
    Verify micro hierarchical P/R/F1 against HiClass when available, and
    against hand-computed values in all cases.

    For the example above: hTP=5, hFP=1, hFN=1
        hP = 5 / (5+1) = 5/6
        hR = 5 / (5+1) = 5/6
        hF1 = 2 * (5/6)^2 / (2 * 5/6) = 5/6
    """

    def _make_scenario(self):
        y_true = ["CL:0000002", "CL:0000002", "CL:0000003"]
        y_pred = ["CL:0000002", "CL:0000002", "CL:0000002"]
        cm = build_confusion_matrix(y_true, y_pred, LABELS)
        return y_true, y_pred, cm

    def test_micro_precision(self):
        _, _, cm = self._make_scenario()
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        np.testing.assert_allclose(m["h_precision_micro"], 5 / 6, rtol=1e-9)

    def test_micro_recall(self):
        _, _, cm = self._make_scenario()
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        np.testing.assert_allclose(m["h_recall_micro"], 5 / 6, rtol=1e-9)

    def test_micro_f1(self):
        _, _, cm = self._make_scenario()
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        np.testing.assert_allclose(m["h_f1_micro"], 5 / 6, rtol=1e-9)

    def test_perfect_micro_is_one(self):
        cm = build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002", "CL:0000003"], LABELS)
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        np.testing.assert_allclose(m["h_f1_micro"], 1.0)

    def test_micro_vs_hiclass(self):
        """Validate hierarchical micro metrics against HiClass."""
        y_true, y_pred, cm = self._make_scenario()
        y_true_paths = [LABEL_PATHS[lbl] for lbl in y_true]
        y_pred_paths = [LABEL_PATHS[lbl] for lbl in y_pred]

        hc_p = hc_metrics.precision(y_true_paths, y_pred_paths)
        hc_r = hc_metrics.recall(y_true_paths, y_pred_paths)
        hc_f1 = hc_metrics.f1(y_true_paths, y_pred_paths)

        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        np.testing.assert_allclose(m["h_precision_micro"], hc_p, rtol=1e-6)
        np.testing.assert_allclose(m["h_recall_micro"], hc_r, rtol=1e-6)
        np.testing.assert_allclose(m["h_f1_micro"], hc_f1, rtol=1e-6)

    def test_micro_vs_hiclass_deeper_hierarchy(self):
        """Validate with grandchild (CL:0000004) predictions using a longer path."""
        y_true = ["CL:0000004", "CL:0000004", "CL:0000002"]
        y_pred = ["CL:0000004", "CL:0000003", "CL:0000001"]
        cm = build_confusion_matrix(y_true, y_pred, LABELS)

        y_true_paths = [LABEL_PATHS[lbl] for lbl in y_true]
        y_pred_paths = [LABEL_PATHS[lbl] for lbl in y_pred]
        hc_p = hc_metrics.precision(y_true_paths, y_pred_paths)
        hc_r = hc_metrics.recall(y_true_paths, y_pred_paths)
        hc_f1 = hc_metrics.f1(y_true_paths, y_pred_paths)

        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        np.testing.assert_allclose(m["h_precision_micro"], hc_p, rtol=1e-6)
        np.testing.assert_allclose(m["h_recall_micro"], hc_r, rtol=1e-6)
        np.testing.assert_allclose(m["h_f1_micro"], hc_f1, rtol=1e-6)


# ---------------------------------------------------------------------------
# compute_hierarchical_f_measure_from_cm — macro F1
# ---------------------------------------------------------------------------


class TestHierarchicalMacro:
    """
    Hand-computed macro F1 (per-node OvR, averaged over nodes with nonzero true support).

    For the base scenario (CL:0000002 x2 correct, CL:0000003 x1 -> CL:0000002):

        Node CL:0000001: tp=3, fp=0, fn=0  -> P=1.0, R=1.0, F1=1.0
        Node CL:0000002: tp=2, fp=1, fn=0  -> P=2/3, R=1.0, F1=4/5=0.8
        Node CL:0000003: tp=0, fp=0, fn=1  -> P=0.0, R=0.0, F1=0.0

        macro F1 = (1.0 + 0.8 + 0.0) / 3 = 0.6
    """

    def test_macro_f1_hand_computed(self):
        cm = build_confusion_matrix(
            ["CL:0000002", "CL:0000002", "CL:0000003"],
            ["CL:0000002", "CL:0000002", "CL:0000002"],
            LABELS,
        )
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        np.testing.assert_allclose(m["h_f1_macro"], 3 / 5, rtol=1e-9)

    def test_macro_f1_perfect_is_one(self):
        cm = build_confusion_matrix(
            ["CL:0000002", "CL:0000003"],
            ["CL:0000002", "CL:0000003"],
            LABELS,
        )
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        np.testing.assert_allclose(m["h_f1_macro"], 1.0)

    def test_macro_uses_only_true_classes(self):
        # Only CL:0000002 appears as a true label; macro is over one class
        cm = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        # Single class with perfect prediction -> macro = 1.0
        np.testing.assert_allclose(m["h_f1_macro"], 1.0)

    def test_macro_symmetric_siblings(self):
        """
        Symmetric case: each sibling misclassified as the other.
            Cell 1: true=CL:0000002, pred=CL:0000003
                T={1,2}, P={1,3}: T∩P={1}->hTP=1, P-T={3}->hFP=1, T-P={2}->hFN=1
            Cell 2: true=CL:0000003, pred=CL:0000002
                T={1,3}, P={1,2}: T∩P={1}->hTP=1, P-T={2}->hFP=1, T-P={3}->hFN=1

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

    def test_zero_division_macro(self):
        cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE, zero_division=0.0)
        np.testing.assert_allclose(m["h_f1_macro"], 0.0)

    def test_h_weighted_keys_present(self):
        cm = build_confusion_matrix(
            ["CL:0000002", "CL:0000002", "CL:0000003"],
            ["CL:0000002", "CL:0000002", "CL:0000002"],
            LABELS,
        )
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        for key in ("h_precision_weighted", "h_recall_weighted", "h_f1_weighted"):
            assert key in m, f"Key {key!r} missing from result"
            assert isinstance(m[key], float), f"{key} should be float"
            assert 0.0 <= m[key] <= 1.0, f"{key}={m[key]} out of [0, 1]"

    def test_h_weighted_perfect(self):
        cm = build_confusion_matrix(
            ["CL:0000002", "CL:0000003"],
            ["CL:0000002", "CL:0000003"],
            LABELS,
        )
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        np.testing.assert_allclose(m["h_precision_weighted"], 1.0)
        np.testing.assert_allclose(m["h_recall_weighted"], 1.0)
        np.testing.assert_allclose(m["h_f1_weighted"], 1.0)

    def test_h_weighted_empty_cm(self):
        cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE, zero_division=0.0)
        np.testing.assert_allclose(m["h_precision_weighted"], 0.0)
        np.testing.assert_allclose(m["h_recall_weighted"], 0.0)
        np.testing.assert_allclose(m["h_f1_weighted"], 0.0)

    def test_h_macro_precision_recall_hand_computed(self):
        """
        Base scenario: CL:0000002 x2 correct, CL:0000003 x1 -> CL:0000002.

        Node CL:0000001: tp=3, fp=0, fn=0 -> P=1.0,  R=1.0
        Node CL:0000002: tp=2, fp=1, fn=0 -> P=2/3,  R=1.0
        Node CL:0000003: tp=0, fp=0, fn=1 -> P=0.0,  R=0.0
        precision_macro = (1.0 + 2/3 + 0.0) / 3 = 5/9
        recall_macro    = (1.0 + 1.0 + 0.0) / 3 = 2/3
        """
        cm = build_confusion_matrix(
            ["CL:0000002", "CL:0000002", "CL:0000003"],
            ["CL:0000002", "CL:0000002", "CL:0000002"],
            LABELS,
        )
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        np.testing.assert_allclose(m["h_precision_macro"], 5 / 9, rtol=1e-9)
        np.testing.assert_allclose(m["h_recall_macro"], 2 / 3, rtol=1e-9)

    def test_h_macro_precision_recall_perfect(self):
        cm = build_confusion_matrix(
            ["CL:0000002", "CL:0000003"],
            ["CL:0000002", "CL:0000003"],
            LABELS,
        )
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE)
        np.testing.assert_allclose(m["h_precision_macro"], 1.0)
        np.testing.assert_allclose(m["h_recall_macro"], 1.0)

    def test_h_macro_precision_recall_empty_cm(self):
        cm = scipy.sparse.csr_matrix((len(LABELS), len(LABELS)), dtype=np.int64)
        m = compute_hierarchical_f_measure_from_cm(cm, LABELS, ONTOLOGY_CACHE, zero_division=0.0)
        np.testing.assert_allclose(m["h_precision_macro"], 0.0)
        np.testing.assert_allclose(m["h_recall_macro"], 0.0)
