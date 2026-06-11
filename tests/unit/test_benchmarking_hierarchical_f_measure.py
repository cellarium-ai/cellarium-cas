import typing as t

import pandas as pd
import pytest

from cellarium.cas.benchmarking.hierarchical_f_measure import (
    compute_hierarchical_f_measure_metrics,
)

# Minimal mock ontology:
#   root (CL:0000000)
#     └── parent (CL:0000001)
#           ├── child_a (CL:0000002)
#           └── child_b (CL:0000003)

MOCK_ONTOLOGY_RESOURCE: t.Dict[str, t.Any] = {
    "cl_names": ["CL:0000000", "CL:0000001", "CL:0000002", "CL:0000003"],
    "cell_ontology_term_id_to_cell_type": {
        "CL:0000000": "root",
        "CL:0000001": "parent",
        "CL:0000002": "child_a",
        "CL:0000003": "child_b",
    },
    "children_dictionary": {
        "CL:0000000": ["CL:0000001"],
        "CL:0000001": ["CL:0000002", "CL:0000003"],
    },
}

ROOT_CHILD_ONTOLOGY_RESOURCE: t.Dict[str, t.Any] = {
    "cl_names": ["CL:0000000", "CL:0000002"],
    "cell_ontology_term_id_to_cell_type": {
        "CL:0000000": "root",
        "CL:0000002": "child_a",
    },
    "children_dictionary": {
        "CL:0000000": ["CL:0000002"],
    },
}


def test_compute_hierarchical_f_measure_length_mismatch_raises():
    predictions = [["CL:0000002"]]

    with pytest.raises(ValueError, match="Length mismatch"):
        compute_hierarchical_f_measure_metrics(
            predictions=predictions,
            ground_truths=["CL:0000002", "CL:0000003"],
            ontology_resource=MOCK_ONTOLOGY_RESOURCE,
        )


def test_compute_hierarchical_f_measure_unrecognized_gt_raises():
    predictions = [["CL:0000002"]]

    with pytest.raises(ValueError, match="not present in the ontology resource"):
        compute_hierarchical_f_measure_metrics(
            predictions=predictions,
            ground_truths=["CL:9999999"],
            ontology_resource=MOCK_ONTOLOGY_RESOURCE,
        )


def test_compute_hierarchical_f_measure_empty_predictions():
    result = compute_hierarchical_f_measure_metrics(
        predictions=[],
        ground_truths=[],
        ontology_resource=MOCK_ONTOLOGY_RESOURCE,
    )

    assert isinstance(result, pd.DataFrame)
    assert result.empty


def test_compute_hierarchical_f_measure_micro_values():
    predictions = [["CL:0000002"], ["CL:0000002"]]

    result = compute_hierarchical_f_measure_metrics(
        predictions=predictions,
        ground_truths=["CL:0000002", "CL:0000003"],
        ontology_resource=MOCK_ONTOLOGY_RESOURCE,
    )

    assert result.loc[0, "top_1_micro_hierarchical_precision"] == pytest.approx(5 / 6)
    assert result.loc[0, "top_1_micro_hierarchical_recall"] == pytest.approx(5 / 6)
    assert result.loc[0, "top_1_micro_hierarchical_f1"] == pytest.approx(5 / 6)


def test_compute_hierarchical_f_measure_top_k_adds_columns_per_k():
    predictions = [["CL:0000003", "CL:0000002"]]

    result = compute_hierarchical_f_measure_metrics(
        predictions=predictions,
        ground_truths=["CL:0000002"],
        ontology_resource=MOCK_ONTOLOGY_RESOURCE,
        top_k=2,
    )

    assert result.loc[0, "top_1_micro_hierarchical_precision"] == pytest.approx(2 / 3)
    assert result.loc[0, "top_1_micro_hierarchical_recall"] == pytest.approx(2 / 3)
    assert result.loc[0, "top_2_micro_hierarchical_precision"] == pytest.approx(3 / 4)
    assert result.loc[0, "top_2_micro_hierarchical_recall"] == pytest.approx(1.0)


def test_compute_hierarchical_f_measure_non_positive_top_k_raises():
    with pytest.raises(ValueError, match="top_k must be positive"):
        compute_hierarchical_f_measure_metrics(
            predictions=[["CL:0000002"]],
            ground_truths=["CL:0000002"],
            ontology_resource=MOCK_ONTOLOGY_RESOURCE,
            top_k=0,
        )


def test_compute_hierarchical_f_measure_summary_columns():
    predictions = [["CL:0000002"]]

    result = compute_hierarchical_f_measure_metrics(
        predictions=predictions,
        ground_truths=["CL:0000002"],
        ontology_resource=MOCK_ONTOLOGY_RESOURCE,
    )

    assert list(result.columns) == [
        "n_cells",
        "top_1_micro_hierarchical_precision",
        "top_1_micro_hierarchical_recall",
        "top_1_micro_hierarchical_f1",
        "top_1_macro_hierarchical_precision",
        "top_1_macro_hierarchical_recall",
        "top_1_macro_hierarchical_f1",
        "top_1_macro_weighted_hierarchical_precision",
        "top_1_macro_weighted_hierarchical_recall",
        "top_1_macro_weighted_hierarchical_f1",
    ]


def test_compute_hierarchical_f_measure_macro_weighted_pools_counts_before_division():
    predictions = [["CL:0000002"], []]

    result = compute_hierarchical_f_measure_metrics(
        predictions=predictions,
        ground_truths=["CL:0000002", "CL:0000002"],
        ontology_resource=ROOT_CHILD_ONTOLOGY_RESOURCE,
        class_level=True,
    )

    row = result.iloc[0]
    assert row["tp"] == pytest.approx(2.0)
    assert row["fp"] == pytest.approx(0.0)
    assert row["fn"] == pytest.approx(2.0)
    assert row["hierarchical_precision"] == pytest.approx(1.0)
    assert row["hierarchical_recall"] == pytest.approx(1 / 2)
    assert row["hierarchical_f1"] == pytest.approx(2 / 3)


def test_compute_hierarchical_f_measure_unrecognized_prediction_raises():
    with pytest.raises(ValueError, match="predicted terms are not present"):
        compute_hierarchical_f_measure_metrics(
            predictions=[["CL:9999999"]],
            ground_truths=["CL:0000002"],
            ontology_resource=MOCK_ONTOLOGY_RESOURCE,
            class_level=True,
        )


def test_compute_hierarchical_f_measure_class_level_columns():
    predictions = [["CL:0000002"]]

    result = compute_hierarchical_f_measure_metrics(
        predictions=predictions,
        ground_truths=["CL:0000002"],
        ontology_resource=MOCK_ONTOLOGY_RESOURCE,
        class_level=True,
    )

    assert list(result.columns) == [
        "ground_truth_class",
        "support",
        "weight",
        "tp",
        "fp",
        "fn",
        "hierarchical_precision",
        "hierarchical_recall",
        "hierarchical_f1",
    ]
