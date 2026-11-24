"""
Unit tests for curated granularity annotation selection.

Synthetic ontology terms:

    CL:0000001  depth 1
    CL:0000002  depth 2
    CL:0000003  depth 3
    CL:0000004  depth 4
    CL:0000005  depth 4 with identical L3/L4 annotations
    CL:0000006  no curated levels
"""

import typing as t

import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from parameterized import parameterized

from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache
from cellarium.cas.postprocessing.cell_ontology.granularity_cache import CellGranularityCache
from cellarium.cas.postprocessing.ontology_aware import (
    CAS_CL_SCORES_ANNDATA_OBSM_KEY,
    CAS_METADATA_ANNDATA_UNS_KEY,
    AggregatedCellOntologyScores,
    CellOntologyScoresAggregationDomain,
    CellOntologyScoresAggregationOp,
    compute_most_granular_top_k_calls_cluster,
    compute_most_granular_top_k_calls_knn,
    compute_most_granular_top_k_calls_single,
    get_most_granular_top_k_calls,
)

TEST_CL_NAMES = [
    "CL:0000001",
    "CL:0000002",
    "CL:0000003",
    "CL:0000004",
    "CL:0000005",
    "CL:0000006",
]
TEST_CL_LABELS = {name: f"label_{name}" for name in TEST_CL_NAMES}

TEST_GRANULARITY_RESOURCE = {
    "CL:0000001": {"l1": "CL:0000001", "l2": None, "l3": None, "l4": None},
    "CL:0000002": {"l1": "CL:0000001", "l2": "CL:0000002", "l3": None, "l4": None},
    "CL:0000003": {"l1": "CL:0000001", "l2": "CL:0000002", "l3": "CL:0000003", "l4": None},
    "CL:0000004": {
        "l1": "CL:0000001",
        "l2": "CL:0000002",
        "l3": "CL:0000003",
        "l4": "CL:0000004",
    },
    "CL:0000005": {
        "l1": "CL:0000001",
        "l2": "CL:0000002",
        "l3": "CL:0000005",
        "l4": "CL:0000005",
    },
    "CL:0000006": {"l1": None, "l2": None, "l3": None, "l4": None},
}


def _make_cl() -> CellOntologyCache:
    resource = {
        "cl_names": TEST_CL_NAMES,
        "cell_ontology_term_id_to_cell_type": TEST_CL_LABELS,
        "children_dictionary": {name: [] for name in TEST_CL_NAMES},
        "shortest_path_lengths_from_cell_root": {name: 0 for name in TEST_CL_NAMES},
        "longest_path_lengths_from_cell_root": {name: 0 for name in TEST_CL_NAMES},
    }
    return CellOntologyCache(resource)


def _aggregated_scores(scores: t.Dict[str, float]) -> AggregatedCellOntologyScores:
    cl_names = list(scores)
    return AggregatedCellOntologyScores(
        aggregation_op=CellOntologyScoresAggregationOp.MEAN,
        aggregation_domain=CellOntologyScoresAggregationDomain.ALL_CELLS,
        threshold=0.0,
        aggregated_scores_c=np.array(list(scores.values())),
        fraction_over_threshold_c=np.ones(len(scores)),
        cl_labels=[TEST_CL_LABELS.get(name, name) for name in cl_names],
        cl_names=cl_names,
    )


def _make_adata(score_matrix: np.ndarray) -> AnnData:
    adata = AnnData(np.zeros((score_matrix.shape[0], 1)))
    adata.obsm[CAS_CL_SCORES_ANNDATA_OBSM_KEY] = sp.csr_matrix(score_matrix)
    adata.uns[CAS_METADATA_ANNDATA_UNS_KEY] = {
        "cl_names": TEST_CL_NAMES,
        "cl_labels": [TEST_CL_LABELS[name] for name in TEST_CL_NAMES],
    }
    return adata


@parameterized.expand(
    [
        ("CL:0000001", 1),
        ("CL:0000002", 2),
        ("CL:0000003", 3),
        ("CL:0000004", 4),
        ("CL:0000005", 4),
        ("CL:0000006", 0),
        ("CL:9999999", 0),
    ]
)
def test_curated_depth(cl_name: str, expected_depth: int):
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    assert cache.get_depth(cl_name) == expected_depth


def test_identical_l3_and_l4_gives_depth_4():
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    assert cache.get_depth("CL:0000005") == 4
    assert cache.get_levels("CL:0000005")["l3"] == "CL:0000005"
    assert cache.get_levels("CL:0000005")["l4"] == "CL:0000005"


def test_missing_levels_and_unknown_values_are_normalized():
    cache = CellGranularityCache(
        {
            "CL:X": {"l1": "CL:A", "l2": "unknown", "l3": None, "l4": ""},
        }
    )
    assert cache.get_depth("CL:X") == 1
    assert cache.get_levels("CL:X") == {"l1": "CL:A", "l2": None, "l3": None, "l4": None}
    assert cache.get_levels("CL:MISSING") == {"l1": None, "l2": None, "l3": None, "l4": None}


def test_deeper_candidates_win_over_higher_scoring_shallow_candidates():
    cl = _make_cl()
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    calls = get_most_granular_top_k_calls(
        aggregated_scores=_aggregated_scores({"CL:0000001": 0.92, "CL:0000003": 0.75, "CL:0000004": 0.45}),
        cl=cl,
        granularity_cache=cache,
        min_acceptable_score=0.2,
        top_k=1,
    )
    assert [call.cl_name for call in calls] == ["CL:0000004"]
    assert calls[0].score == 0.45
    assert calls[0].levels == TEST_GRANULARITY_RESOURCE["CL:0000004"]


def test_highest_score_wins_within_same_curated_depth():
    cl = _make_cl()
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    calls = get_most_granular_top_k_calls(
        aggregated_scores=_aggregated_scores({"CL:0000004": 0.4, "CL:0000005": 0.8}),
        cl=cl,
        granularity_cache=cache,
        min_acceptable_score=0.2,
        top_k=2,
    )
    assert [call.cl_name for call in calls] == ["CL:0000005", "CL:0000004"]


def test_fallback_when_no_l4_candidate_meets_threshold():
    cl = _make_cl()
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    calls = get_most_granular_top_k_calls(
        aggregated_scores=_aggregated_scores({"CL:0000003": 0.75, "CL:0000004": 0.1}),
        cl=cl,
        granularity_cache=cache,
        min_acceptable_score=0.2,
        top_k=1,
    )
    assert [call.cl_name for call in calls] == ["CL:0000003"]


def test_terms_absent_from_lookup_or_with_depth_zero_are_excluded():
    cl = _make_cl()
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    calls = get_most_granular_top_k_calls(
        aggregated_scores=_aggregated_scores({"CL:0000006": 0.99}),
        cl=cl,
        granularity_cache=cache,
        min_acceptable_score=0.2,
        top_k=3,
    )
    assert calls == []

    calls = get_most_granular_top_k_calls(
        aggregated_scores=_aggregated_scores({"CL:0000002": 0.99, "CL:NOT_IN_LOOKUP": 0.99}),
        cl=cl,
        granularity_cache=cache,
        min_acceptable_score=0.2,
        top_k=3,
    )
    assert [call.cl_name for call in calls] == ["CL:0000002"]


def test_deterministic_tie_break_by_cl_name():
    cl = _make_cl()
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    calls = get_most_granular_top_k_calls(
        aggregated_scores=_aggregated_scores({"CL:0000004": 0.5, "CL:0000005": 0.5}),
        cl=cl,
        granularity_cache=cache,
        min_acceptable_score=0.2,
        top_k=2,
    )
    assert [call.cl_name for call in calls] == ["CL:0000004", "CL:0000005"]


def test_fill_missing_levels_flag_controls_undefined_levels():
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    cl = _make_cl()
    for fill in (True, False):
        adata = _make_adata(np.array([[0.92, 0.75, 0.1, 0.05, 0.05, 0.0]]))
        compute_most_granular_top_k_calls_single(
            adata=adata,
            cl=cl,
            granularity_cache=cache,
            min_acceptable_score=0.8,
            top_k=1,
            fill_missing_levels=fill,
        )
        assert adata.obs["cas_cell_type_l1_name_1"].iloc[0] == "CL:0000001"
        assert (adata.obs["cas_cell_type_l4_name_1"].iloc[0] == "CL:0000001") is fill
        assert pd.isna(adata.obs["cas_cell_type_l4_label_1"].iloc[0]) is (not fill)


def test_single_cell_outputs_and_missing_call_padding():
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    cl = _make_cl()
    # Only the depth-2 term clears the threshold, so rank 1 is filled and rank 2 is missing.
    scores = np.array([[0.92, 0.75, 0.1, 0.05, 0.05, 0.0]])
    adata = _make_adata(scores)

    compute_most_granular_top_k_calls_single(
        adata=adata,
        cl=cl,
        granularity_cache=cache,
        min_acceptable_score=0.8,
        top_k=2,
        obs_prefix="cas_cell_type",
    )

    assert adata.obs["cas_cell_type_name_1"].iloc[0] == "CL:0000001"
    assert adata.obs["cas_cell_type_score_1"].iloc[0] == 0.92
    assert adata.obs["cas_cell_type_l1_name_1"].iloc[0] == "CL:0000001"
    assert adata.obs["cas_cell_type_l1_label_1"].iloc[0] == cl.cl_names_to_labels_map["CL:0000001"]
    # CL:0000001 only defines L1, so L2-L4 inherit it.
    assert adata.obs["cas_cell_type_l4_name_1"].iloc[0] == "CL:0000001"
    assert adata.obs["cas_cell_type_l4_label_1"].iloc[0] == cl.cl_names_to_labels_map["CL:0000001"]
    assert pd.isna(adata.obs["cas_cell_type_name_2"].iloc[0])
    assert pd.isna(adata.obs["cas_cell_type_score_2"].iloc[0])
    assert pd.isna(adata.obs["cas_cell_type_l4_name_2"].iloc[0])


def test_no_eligible_candidates_do_not_fabricate_confidence():
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    cl = _make_cl()
    adata = _make_adata(np.array([[0.1, 0.1, 0.1, 0.1, 0.1, 0.1]]))

    compute_most_granular_top_k_calls_single(
        adata=adata,
        cl=cl,
        granularity_cache=cache,
        min_acceptable_score=0.2,
        top_k=3,
    )

    for rank in range(1, 4):
        assert pd.isna(adata.obs[f"cas_cell_type_name_{rank}"].iloc[0])
        assert pd.isna(adata.obs[f"cas_cell_type_score_{rank}"].iloc[0])
        assert pd.isna(adata.obs[f"cas_cell_type_l1_name_{rank}"].iloc[0])


def test_cluster_outputs_use_aggregated_scores():
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    cl = _make_cl()

    # Cluster A has mean depth-2 score 0.5; cluster B has mean 0.3.
    scores = np.zeros((4, len(TEST_CL_NAMES)))
    scores[:, 1] = [0.4, 0.6, 0.2, 0.4]
    adata = _make_adata(scores)
    adata.obs["cluster_label"] = pd.Categorical(["A", "A", "B", "B"])

    compute_most_granular_top_k_calls_cluster(
        adata=adata,
        cl=cl,
        granularity_cache=cache,
        min_acceptable_score=0.2,
        cluster_label_obs_column="cluster_label",
        top_k=2,
        obs_prefix="cas_cluster_type",
    )

    assert adata.obs["cas_cluster_type_name_1"].tolist() == ["CL:0000002"] * 4
    assert np.allclose(adata.obs["cas_cluster_type_score_1"].to_numpy(), [0.5, 0.5, 0.3, 0.3])
    assert adata.obs["cas_cluster_type_l2_name_1"].tolist() == ["CL:0000002"] * 4
    assert pd.isna(adata.obs["cas_cluster_type_name_2"].iloc[0])


def test_knn_outputs_use_neighborhood_aggregation():
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    cl = _make_cl()

    scores = np.zeros((3, len(TEST_CL_NAMES)))
    scores[:, 2] = [0.3, 0.6, 0.9]
    adata = _make_adata(scores)

    # Each cell is connected to the other two cells.
    rows = np.repeat(np.arange(3), 2)
    cols = np.array([1, 2, 0, 2, 0, 1])
    adata.obsp["distances"] = sp.csr_matrix((np.ones(6), (rows, cols)), shape=(3, 3))

    compute_most_granular_top_k_calls_knn(
        adata=adata,
        cl=cl,
        granularity_cache=cache,
        min_acceptable_score=0.2,
        k_neighbors=2,
        top_k=1,
        obs_prefix="cas_knn_cell_type",
    )

    expected_score = np.mean(scores[:, 2])
    assert adata.obs["cas_knn_cell_type_name_1"].tolist() == ["CL:0000003"] * 3
    assert np.allclose(pd.to_numeric(adata.obs["cas_knn_cell_type_score_1"]).to_numpy(), expected_score)
    assert adata.obs["cas_knn_cell_type_l3_name_1"].tolist() == ["CL:0000003"] * 3
