"""
Unit tests for build_confusion_matrix, aggregate_confusion_matrices, and save/load round-trip.
"""

import numpy as np
import pytest
import scipy.sparse

from cellarium.cas.benchmarking.confusion_matrix import (
    aggregate_confusion_matrices,
    build_confusion_matrix,
    load_confusion_matrix,
    save_confusion_matrix,
)
from tests.unit.benchmarking._fixtures import LABELS, MOCK_ONTOLOGY_RESOURCE

# ---------------------------------------------------------------------------
# build_confusion_matrix
# ---------------------------------------------------------------------------


def test_build_cm_returns_sparse_matrix():
    cm = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
    assert scipy.sparse.issparse(cm)


def test_build_cm_shape_equals_label_count():
    cm = build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002", "CL:0000002"], LABELS)
    assert cm.shape == (len(LABELS), len(LABELS))


def test_build_cm_correct_prediction_on_diagonal():
    # CL:0000002 (index 1) predicted correctly
    cm = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
    dense = cm.toarray()
    assert dense[1, 1] == 1  # row=true=CL:0000002, col=pred=CL:0000002


def test_build_cm_off_diagonal_for_wrong_prediction():
    # true=CL:0000003 (idx 2), pred=CL:0000002 (idx 1)
    cm = build_confusion_matrix(["CL:0000003"], ["CL:0000002"], LABELS)
    dense = cm.toarray()
    assert dense[2, 1] == 1
    assert dense[1, 1] == 0


def test_build_cm_multiple_cells_correct_counts():
    y_true = ["CL:0000002", "CL:0000002", "CL:0000003"]
    y_pred = ["CL:0000002", "CL:0000002", "CL:0000002"]
    cm = build_confusion_matrix(y_true, y_pred, LABELS)
    dense = cm.toarray()
    assert dense[1, 1] == 2  # 2 correct for CL:0000002
    assert dense[2, 1] == 1  # 1 wrong: CL:0000003 -> CL:0000002


def test_build_cm_raises_on_unmapped_gt_label():
    with pytest.raises(ValueError, match="not in the ontology"):
        build_confusion_matrix(["CL:9999999"], ["CL:0000002"], LABELS)


def test_build_cm_raises_on_unmapped_pred_label():
    with pytest.raises(ValueError, match="not in the ontology"):
        build_confusion_matrix(["CL:0000002"], ["CL:9999999"], LABELS)


def test_build_cm_raises_on_length_mismatch():
    with pytest.raises(ValueError, match="same length"):
        build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002"], LABELS)


def test_build_cm_full_universe_produces_correct_shape():
    full_labels = MOCK_ONTOLOGY_RESOURCE["cl_names"]
    cm = build_confusion_matrix(["CL:0000001"], ["CL:0000001"], full_labels)
    assert cm.shape == (len(full_labels), len(full_labels))


# ---------------------------------------------------------------------------
# aggregate_confusion_matrices
# ---------------------------------------------------------------------------


def test_aggregate_cm_sums_two_matrices():
    cm1 = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
    cm2 = build_confusion_matrix(["CL:0000003"], ["CL:0000002"], LABELS)
    agg = aggregate_confusion_matrices([cm1, cm2])
    dense = agg.toarray()
    assert dense[1, 1] == 1  # one correct CL:0000002
    assert dense[2, 1] == 1  # one wrong: CL:0000003 -> CL:0000002


def test_aggregate_cm_preserves_shape():
    cm1 = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
    cm2 = build_confusion_matrix(["CL:0000001"], ["CL:0000001"], LABELS)
    agg = aggregate_confusion_matrices([cm1, cm2])
    assert agg.shape == cm1.shape


def test_aggregate_cm_single_matrix_unchanged():
    cm = build_confusion_matrix(["CL:0000002", "CL:0000003"], ["CL:0000002", "CL:0000003"], LABELS)
    agg = aggregate_confusion_matrices([cm])
    np.testing.assert_array_equal(agg.toarray(), cm.toarray())


def test_aggregate_cm_raises_on_empty_list():
    with pytest.raises(ValueError, match="empty"):
        aggregate_confusion_matrices([])


def test_aggregate_cm_raises_on_shape_mismatch():
    cm1 = build_confusion_matrix(["CL:0000002"], ["CL:0000002"], LABELS)
    cm2 = scipy.sparse.csr_matrix(np.zeros((3, 3), dtype=np.int64))
    with pytest.raises(ValueError, match="shape"):
        aggregate_confusion_matrices([cm1, cm2])


# ---------------------------------------------------------------------------
# save / load confusion matrix (round-trip)
# ---------------------------------------------------------------------------


def test_save_load_cm_roundtrip(tmp_path):
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


def test_load_cm_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="confusion_matrix.npz"):
        load_confusion_matrix(tmp_path)
