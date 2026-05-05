import typing as t

import anndata
import numpy as np
import pandas as pd
import pytest

from cellarium.cas.benchmarking.flat import (
    _effective_predictions_at_k,
    compute_flat_metrics,
    extract_predictions_from_adata,
)

# Fixtures


@pytest.fixture
def perfect_predictions() -> t.Tuple[t.List[str], t.List[t.List[str]]]:
    """All top-1 predictions are correct."""
    gts = ["A", "B", "C", "A"]
    preds = [["A", "B", "C"], ["B", "A", "C"], ["C", "A", "B"], ["A", "C", "B"]]
    return gts, preds


@pytest.fixture
def all_wrong_predictions() -> t.Tuple[t.List[str], t.List[t.List[str]]]:
    """Top-1 predictions are always wrong; GT never appears in any rank."""
    gts = ["A", "B", "C"]
    preds = [["B", "C"], ["A", "C"], ["A", "B"]]
    return gts, preds


@pytest.fixture
def top2_hit_predictions() -> t.Tuple[t.List[str], t.List[t.List[str]]]:
    """GT appears at rank-2 (not rank-1) for all cells."""
    gts = ["A", "B", "C"]
    preds = [["B", "A"], ["C", "B"], ["A", "C"]]
    return gts, preds


@pytest.fixture
def sample_adata() -> anndata.AnnData:
    """AnnData with CAS prediction columns already inserted."""
    obs = pd.DataFrame(
        {
            "cas_cell_type_label_1": ["T cell", "B cell", "NK cell"],
            "cas_cell_type_label_2": ["NK cell", "T cell", "B cell"],
            "cas_cell_type_label_3": ["B cell", "NK cell", "T cell"],
        },
        index=["cell_0", "cell_1", "cell_2"],
    )
    return anndata.AnnData(X=np.zeros((3, 5)), obs=obs)


# Tests: extract_predictions_from_adata


def test_extract_predictions_returns_correct_shape(sample_adata):
    result = extract_predictions_from_adata(sample_adata, column_prefix="cas_cell_type_label", top_k=3)
    assert len(result) == 3
    assert all(len(cell_preds) == 3 for cell_preds in result)


def test_extract_predictions_correct_values(sample_adata):
    result = extract_predictions_from_adata(sample_adata, column_prefix="cas_cell_type_label", top_k=2)
    assert result[0] == ["T cell", "NK cell"]
    assert result[1] == ["B cell", "T cell"]
    assert result[2] == ["NK cell", "B cell"]


def test_extract_predictions_top_k_capped(sample_adata):
    result = extract_predictions_from_adata(sample_adata, column_prefix="cas_cell_type_label", top_k=1)
    assert result == [["T cell"], ["B cell"], ["NK cell"]]


def test_extract_predictions_missing_columns_truncates(sample_adata):
    # Only top_k=2 columns exist; requesting top_k=5 should return what's available
    result = extract_predictions_from_adata(sample_adata, column_prefix="cas_cell_type_label", top_k=5)
    assert all(len(cell_preds) == 3 for cell_preds in result)


# Tests: _effective_predictions_at_k


def test_effective_preds_gt_in_top_k():
    gts = ["A", "B"]
    preds = [["X", "A"], ["Y", "B"]]
    _, effective = _effective_predictions_at_k(gts, preds, k=2)
    assert effective == ["A", "B"]


def test_effective_preds_gt_not_in_top_k_falls_back_to_top1():
    gts = ["A", "B"]
    preds = [["X", "Y"], ["Z", "W"]]
    _, effective = _effective_predictions_at_k(gts, preds, k=2)
    assert effective == ["X", "Z"]


def test_effective_preds_empty_predictions():
    gts = ["A"]
    preds = [[]]
    _, effective = _effective_predictions_at_k(gts, preds, k=1)
    assert effective == [""]


# Tests: compute_flat_metrics — summary mode


def test_compute_flat_metrics_perfect_predictions_accuracy_1(perfect_predictions):
    gts, preds = perfect_predictions
    df = compute_flat_metrics(gts, preds, top_k=1)
    assert len(df) == 1
    assert df.loc[0, "top_1_accuracy"] == pytest.approx(1.0)


def test_compute_flat_metrics_all_wrong_accuracy_0_at_k1(all_wrong_predictions):
    gts, preds = all_wrong_predictions
    df = compute_flat_metrics(gts, preds, top_k=1)
    assert df.loc[0, "top_1_accuracy"] == pytest.approx(0.0)


def test_compute_flat_metrics_top2_hit_accuracy_at_k2(top2_hit_predictions):
    gts, preds = top2_hit_predictions
    df = compute_flat_metrics(gts, preds, top_k=2)
    # At k=1: all wrong (GT is at rank-2)
    assert df.loc[0, "top_1_accuracy"] == pytest.approx(0.0)
    # At k=2: all correct
    assert df.loc[0, "top_2_accuracy"] == pytest.approx(1.0)


def test_compute_flat_metrics_returns_all_k_columns(perfect_predictions):
    gts, preds = perfect_predictions
    df = compute_flat_metrics(gts, preds, top_k=3)
    assert len(df) == 1
    for k in range(1, 4):
        assert f"top_{k}_accuracy" in df.columns


def test_compute_flat_metrics_columns_present(perfect_predictions):
    gts, preds = perfect_predictions
    df = compute_flat_metrics(gts, preds, top_k=1)
    expected_cols = {
        "n_cells",
        "top_1_accuracy",
        "top_1_macro_f1",
        "top_1_weighted_f1",
        "top_1_macro_precision",
        "top_1_weighted_precision",
        "top_1_macro_recall",
        "top_1_weighted_recall",
    }
    assert expected_cols.issubset(set(df.columns))


def test_compute_flat_metrics_top_k_inferred_when_none(perfect_predictions):
    gts, preds = perfect_predictions
    df = compute_flat_metrics(gts, preds)  # top_k=None → infer from min pred length
    assert len(df) == 1
    # all pred lists have length 3, so top_1 through top_3 columns should exist
    for k in range(1, 4):
        assert f"top_{k}_accuracy" in df.columns


def test_compute_flat_metrics_length_mismatch_raises():
    with pytest.raises(ValueError, match="Length mismatch"):
        compute_flat_metrics(["A", "B"], [["A"]])


def test_compute_flat_metrics_empty_returns_empty_df():
    df = compute_flat_metrics([], [])
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 0


# Tests: compute_flat_metrics — cell-level mode


def test_compute_flat_metrics_cell_level_shape(perfect_predictions):
    gts, preds = perfect_predictions
    df = compute_flat_metrics(gts, preds, top_k=2, cell_level=True)
    # one row per cell
    assert len(df) == 4


def test_compute_flat_metrics_cell_level_columns(perfect_predictions):
    gts, preds = perfect_predictions
    df = compute_flat_metrics(gts, preds, top_k=1, cell_level=True)
    assert set(df.columns) == {"cell_index", "ground_truth", "top_1_effective_prediction", "top_1_correct"}


def test_compute_flat_metrics_cell_level_correct_flags(perfect_predictions):
    gts, preds = perfect_predictions
    df = compute_flat_metrics(gts, preds, top_k=1, cell_level=True)
    # All top-1 are correct for perfect_predictions
    assert df["top_1_correct"].all()
