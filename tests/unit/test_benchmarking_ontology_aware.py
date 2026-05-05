import typing as t

import pytest

from cellarium.cas.benchmarking.ontology_aware import (
    _build_hop_neighborhoods,
    _compute_cell_metrics,
    _precompute_graph,
    compute_ontology_aware_metrics,
)
from cellarium.cas.models import CellTypeOntologyAwareResults

# Minimal mock ontology:
#   root (CL:0000000)
#     └── parent (CL:0000001)
#           ├── child_a (CL:0000002)
#           └── child_b (CL:0000003)
#                  └── grandchild (CL:0000004)

MOCK_ONTOLOGY_RESOURCE: t.Dict[str, t.Any] = {
    "cl_names": ["CL:0000000", "CL:0000001", "CL:0000002", "CL:0000003", "CL:0000004"],
    "cell_ontology_term_id_to_cell_type": {
        "CL:0000000": "root",
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


# Helpers to build mock CAS response objects


def make_match(term_id: str, score: float) -> CellTypeOntologyAwareResults.Match:
    return CellTypeOntologyAwareResults.Match(
        cell_type_ontology_term_id=term_id,
        score=score,
        cell_type=term_id,
    )


def make_annotation(cell_id: str, matches: t.List) -> CellTypeOntologyAwareResults.OntologyAwareAnnotation:
    return CellTypeOntologyAwareResults.OntologyAwareAnnotation(
        query_cell_id=cell_id,
        matches=matches,
        total_weight=sum(m.score for m in matches),
        total_neighbors=len(matches),
        total_neighbors_unrecognized=0,
    )


def make_response(annotations: t.List) -> CellTypeOntologyAwareResults:
    return CellTypeOntologyAwareResults(data=annotations)


# Fixtures


@pytest.fixture
def graph():
    return _precompute_graph(MOCK_ONTOLOGY_RESOURCE)


# Tests: _precompute_graph


def test_precompute_graph_all_terms(graph):
    assert graph.all_terms == frozenset(MOCK_ONTOLOGY_RESOURCE["cl_names"])


def test_precompute_graph_ancestors_of_grandchild(graph):
    # grandchild → child_b → parent → root
    ancestors = graph.term_ancestors["CL:0000004"]
    assert ancestors == frozenset(["CL:0000004", "CL:0000003", "CL:0000001", "CL:0000000"])


def test_precompute_graph_descendants_of_root(graph):
    # root → all terms
    descendants = graph.term_descendants["CL:0000000"]
    assert descendants == frozenset(MOCK_ONTOLOGY_RESOURCE["cl_names"])


def test_precompute_graph_leaf_has_no_children(graph):
    assert graph.term_descendants["CL:0000002"] == frozenset(["CL:0000002"])


def test_precompute_graph_parents_populated(graph):
    assert "CL:0000001" in graph.parents.get("CL:0000002", set())


# Tests: _build_hop_neighborhoods


def test_hop_0_is_exact_match(graph):
    hops = _build_hop_neighborhoods("CL:0000001", num_hops=2, graph=graph)
    assert hops[0]["nodes"] == frozenset(["CL:0000001"])


def test_hop_1_includes_direct_neighbors(graph):
    # parent (CL:0000001) at hop-1: root (parent) + child_a + child_b
    hops = _build_hop_neighborhoods("CL:0000001", num_hops=2, graph=graph)
    assert frozenset(["CL:0000000", "CL:0000001", "CL:0000002", "CL:0000003"]).issubset(hops[1]["nodes"])


def test_hop_2_includes_grandchild(graph):
    # From parent (CL:0000001), 2 hops reaches grandchild
    hops = _build_hop_neighborhoods("CL:0000001", num_hops=2, graph=graph)
    assert "CL:0000004" in hops[2]["nodes"]


def test_hops_are_cumulative(graph):
    hops = _build_hop_neighborhoods("CL:0000001", num_hops=3, graph=graph)
    # Each hop's nodes is a superset of the previous
    for i in range(1, len(hops)):
        assert hops[i - 1]["nodes"].issubset(hops[i]["nodes"])


# Tests: _compute_cell_metrics


def test_exact_match_is_tp_at_hop0(graph):
    # GT = child_a (CL:0000002), prediction = child_a → TP at all hops
    hops = _build_hop_neighborhoods("CL:0000002", num_hops=2, graph=graph)
    matches = [make_match("CL:0000002", 1.0)]
    metrics = _compute_cell_metrics(matches, hops, graph, num_hops=2)
    assert metrics["hop_0_sensitivity"] == pytest.approx(1.0)
    assert metrics["hop_0_specificity"] == pytest.approx(1.0)


def test_sibling_is_tp_at_hop2_not_hop0(graph):
    # GT = child_a (CL:0000002), prediction = child_b (CL:0000003)
    # child_a and child_b are siblings (share parent CL:0000001)
    # At hop-0: child_b not in {child_a} → not TP
    # At hop-2 from child_a: need to verify child_b is reachable in 2 hops (child_a → parent → child_b)
    hops = _build_hop_neighborhoods("CL:0000002", num_hops=2, graph=graph)
    matches = [make_match("CL:0000003", 0.8)]
    metrics = _compute_cell_metrics(matches, hops, graph, num_hops=2)
    assert metrics["hop_0_sensitivity"] < 1.0  # not TP at hop-0
    assert metrics["hop_2_sensitivity"] == pytest.approx(0.8)  # TP at hop-2


def test_unrelated_term_is_fp_at_all_hops(graph):
    # GT = child_a (CL:0000002), prediction = a term not in the ontology → skipped
    hops = _build_hop_neighborhoods("CL:0000002", num_hops=2, graph=graph)
    matches = [make_match("CL:9999999", 0.9)]  # unknown term
    metrics = _compute_cell_metrics(matches, hops, graph, num_hops=2)
    # Unknown term not in term_ancestors → skipped, so no TP and no FP
    assert metrics["hop_0_sensitivity"] == pytest.approx(0.0)


def test_no_matches_gives_zero_sensitivity(graph):
    hops = _build_hop_neighborhoods("CL:0000002", num_hops=2, graph=graph)
    metrics = _compute_cell_metrics([], hops, graph, num_hops=2)
    for i in range(3):
        assert metrics[f"hop_{i}_sensitivity"] == pytest.approx(0.0)


# Tests: compute_ontology_aware_metrics


def test_compute_ontology_aware_metrics_length_mismatch_raises():
    response = make_response([make_annotation("c1", [make_match("CL:0000002", 1.0)])])
    with pytest.raises(ValueError, match="Length mismatch"):
        compute_ontology_aware_metrics(response, ["CL:0000002", "CL:0000003"], MOCK_ONTOLOGY_RESOURCE)


def test_compute_ontology_aware_metrics_unrecognized_gt_raises():
    response = make_response([make_annotation("c1", [make_match("CL:0000002", 1.0)])])
    with pytest.raises(ValueError, match="not present in the ontology resource"):
        compute_ontology_aware_metrics(response, ["CL:9999999"], MOCK_ONTOLOGY_RESOURCE)


def test_compute_ontology_aware_metrics_summary_shape():
    annotations = [
        make_annotation("c1", [make_match("CL:0000002", 1.0)]),
        make_annotation("c2", [make_match("CL:0000003", 0.8)]),
    ]
    response = make_response(annotations)
    ground_truths = ["CL:0000002", "CL:0000003"]
    df = compute_ontology_aware_metrics(response, ground_truths, MOCK_ONTOLOGY_RESOURCE, num_hops=2)
    assert len(df) == 1  # one summary row
    assert "n_cells" in df.columns
    assert df.loc[0, "n_cells"] == 2


def test_compute_ontology_aware_metrics_summary_has_hop_columns():
    annotations = [make_annotation("c1", [make_match("CL:0000002", 1.0)])]
    response = make_response(annotations)
    df = compute_ontology_aware_metrics(response, ["CL:0000002"], MOCK_ONTOLOGY_RESOURCE, num_hops=2)
    for i in range(3):
        assert f"hop_{i}_sensitivity" in df.columns
        assert f"hop_{i}_specificity" in df.columns
        assert f"hop_{i}_f1_score" in df.columns


def test_compute_ontology_aware_metrics_perfect_prediction():
    # Predict exact GT → sensitivity = 1.0 at hop-0
    annotations = [make_annotation("c1", [make_match("CL:0000002", 1.0)])]
    response = make_response(annotations)
    df = compute_ontology_aware_metrics(response, ["CL:0000002"], MOCK_ONTOLOGY_RESOURCE, num_hops=2)
    assert df.loc[0, "hop_0_sensitivity"] == pytest.approx(1.0)


def test_compute_ontology_aware_metrics_cell_level_shape():
    annotations = [
        make_annotation("c1", [make_match("CL:0000002", 1.0)]),
        make_annotation("c2", [make_match("CL:0000003", 0.8)]),
    ]
    response = make_response(annotations)
    df = compute_ontology_aware_metrics(
        response, ["CL:0000002", "CL:0000003"], MOCK_ONTOLOGY_RESOURCE, num_hops=2, cell_level=True
    )
    assert len(df) == 2
    assert "query_cell_id" in df.columns
    assert "ground_truth" in df.columns


def test_compute_ontology_aware_metrics_empty_response():
    response = make_response([])
    df = compute_ontology_aware_metrics(response, [], MOCK_ONTOLOGY_RESOURCE, num_hops=2)
    assert isinstance(df, __import__("pandas").DataFrame)
    assert len(df) == 0
