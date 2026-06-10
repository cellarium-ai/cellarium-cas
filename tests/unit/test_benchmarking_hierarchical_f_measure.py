import typing as t

import pandas as pd
import pytest

from cellarium.cas.benchmarking.hierarchical_f_measure import compute_hierarchical_f_measure_metrics
from cellarium.cas.models import CellTypeOntologyAwareResults

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


def make_match(term_id: str, score: float = 1.0) -> CellTypeOntologyAwareResults.Match:
    return CellTypeOntologyAwareResults.Match(
        cell_type_ontology_term_id=term_id,
        score=score,
        cell_type=term_id,
    )


def make_annotation(cell_id: str, matches: t.List[CellTypeOntologyAwareResults.Match]):
    return CellTypeOntologyAwareResults.OntologyAwareAnnotation(
        query_cell_id=cell_id,
        matches=matches,
        total_weight=sum(m.score for m in matches),
        total_neighbors=len(matches),
        total_neighbors_unrecognized=0,
    )


def make_response(annotations: t.List[CellTypeOntologyAwareResults.OntologyAwareAnnotation]):
    return CellTypeOntologyAwareResults(data=annotations)


def test_compute_hierarchical_f_measure_length_mismatch_raises():
    response = make_response([make_annotation("c1", [make_match("CL:0000002")])])

    with pytest.raises(ValueError, match="Length mismatch"):
        compute_hierarchical_f_measure_metrics(
            response=response,
            ground_truths=["CL:0000002", "CL:0000003"],
            ontology_resource=MOCK_ONTOLOGY_RESOURCE,
        )


def test_compute_hierarchical_f_measure_unrecognized_gt_raises():
    response = make_response([make_annotation("c1", [make_match("CL:0000002")])])

    with pytest.raises(ValueError, match="not present in the ontology resource"):
        compute_hierarchical_f_measure_metrics(
            response=response,
            ground_truths=["CL:9999999"],
            ontology_resource=MOCK_ONTOLOGY_RESOURCE,
        )


def test_compute_hierarchical_f_measure_empty_response():
    response = make_response([])

    result = compute_hierarchical_f_measure_metrics(
        response=response,
        ground_truths=[],
        ontology_resource=MOCK_ONTOLOGY_RESOURCE,
    )

    assert isinstance(result, pd.DataFrame)
    assert result.empty


def test_compute_hierarchical_f_measure_micro_values():
    response = make_response(
        [
            make_annotation("c1", [make_match("CL:0000000"), make_match("CL:0000001"), make_match("CL:0000002")]),
            make_annotation("c2", [make_match("CL:0000000"), make_match("CL:0000001"), make_match("CL:0000002")]),
        ]
    )

    result = compute_hierarchical_f_measure_metrics(
        response=response,
        ground_truths=["CL:0000002", "CL:0000003"],
        ontology_resource=MOCK_ONTOLOGY_RESOURCE,
    )

    assert result.loc[0, "hierarchical_precision"] == pytest.approx(5 / 6)
    assert result.loc[0, "hierarchical_recall"] == pytest.approx(5 / 6)
    assert result.loc[0, "hierarchical_f1"] == pytest.approx(5 / 6)


def test_compute_hierarchical_f_measure_macro_weighted_pools_counts_before_division():
    response = make_response(
        [
            make_annotation("c1", [make_match("CL:0000000"), make_match("CL:0000002")]),
            make_annotation("c2", [make_match("unknown")]),
        ]
    )

    result = compute_hierarchical_f_measure_metrics(
        response=response,
        ground_truths=["CL:0000002", "CL:0000002"],
        ontology_resource=ROOT_CHILD_ONTOLOGY_RESOURCE,
        class_level=True,
    )

    row = result.iloc[0]
    assert row["tp"] == 2.0
    assert row["fp"] == 1.0
    assert row["fn"] == 2.0
    assert row["hierarchical_precision"] == pytest.approx(2 / 3)
    assert row["hierarchical_recall"] == pytest.approx(1 / 2)
    assert row["hierarchical_f1"] == pytest.approx(4 / 7)


def test_compute_hierarchical_f_measure_unknown_prediction_counts_as_fp():
    response = make_response([make_annotation("c1", [make_match("CL:9999999")])])

    result = compute_hierarchical_f_measure_metrics(
        response=response,
        ground_truths=["CL:0000002"],
        ontology_resource=MOCK_ONTOLOGY_RESOURCE,
        class_level=True,
    )

    assert result.loc[0, "fp"] == 1.0


def test_compute_hierarchical_f_measure_class_level_columns():
    response = make_response([make_annotation("c1", [make_match("CL:0000002")])])

    result = compute_hierarchical_f_measure_metrics(
        response=response,
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
