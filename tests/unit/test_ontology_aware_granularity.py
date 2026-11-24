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
    MISSING_GRANULARITY_CELL_TYPE_CALL,
    AggregatedCellOntologyScores,
    CellOntologyScoresAggregationDomain,
    CellOntologyScoresAggregationOp,
    compute_most_granular_top_k_calls_cluster,
    compute_most_granular_top_k_calls_knn,
    compute_most_granular_top_k_calls_single,
    get_most_granular_call,
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


def _select(scores: t.Dict[str, float], cl: t.Optional[CellOntologyCache] = None, min_score: float = 0.2):
    return get_most_granular_call(
        aggregated_scores=_aggregated_scores(scores),
        cl=cl or _make_cl(),
        granularity_cache=CellGranularityCache(TEST_GRANULARITY_RESOURCE),
        min_acceptable_score=min_score,
    )


def test_deeper_candidate_wins_over_higher_scoring_shallow_candidate():
    call = _select({"CL:0000001": 0.92, "CL:0000003": 0.75, "CL:0000004": 0.45})
    assert call.cl_name == "CL:0000004"
    assert call.score == 0.45
    assert call.levels == TEST_GRANULARITY_RESOURCE["CL:0000004"]


def test_highest_score_wins_within_same_curated_depth():
    assert _select({"CL:0000004": 0.4, "CL:0000005": 0.8}).cl_name == "CL:0000005"


def test_fallback_when_no_l4_candidate_meets_threshold():
    assert _select({"CL:0000003": 0.75, "CL:0000004": 0.1}).cl_name == "CL:0000003"


def test_terms_absent_from_lookup_or_with_depth_zero_are_excluded():
    assert _select({"CL:0000006": 0.99}) is MISSING_GRANULARITY_CELL_TYPE_CALL
    assert _select({"CL:0000002": 0.99, "CL:NOT_IN_LOOKUP": 0.99}).cl_name == "CL:0000002"


def test_no_eligible_candidate_returns_missing_call():
    call = _select({"CL:0000004": 0.1})
    assert call is MISSING_GRANULARITY_CELL_TYPE_CALL
    assert call.score is None and call.cl_name is None and call.level_scores == {}


def test_deterministic_tie_break_by_cl_name():
    assert _select({"CL:0000005": 0.5, "CL:0000004": 0.5}).cl_name == "CL:0000004"


def test_level_scores_hold_aggregated_scores_of_level_terms_even_below_threshold():
    call = _select({"CL:0000001": 0.1, "CL:0000002": 0.3, "CL:0000003": 0.6, "CL:0000004": 0.45, "CL:0000005": 0.9})
    assert call.cl_name == "CL:0000005"
    assert call.level_scores == {"CL:0000001": 0.1, "CL:0000002": 0.3, "CL:0000005": 0.9}


def _run_single(scores_row, min_score: float = 0.2, **kwargs) -> AnnData:
    adata = _make_adata(np.array([scores_row]))
    compute_most_granular_top_k_calls_single(
        adata=adata,
        cl=_make_cl(),
        granularity_cache=CellGranularityCache(TEST_GRANULARITY_RESOURCE),
        min_acceptable_score=min_score,
        **kwargs,
    )
    return adata


def test_single_positions_run_from_most_to_least_granular_level():
    # CL:0000004 is selected (depth 4): position 1 = l4 ... position 4 = l1.
    adata = _run_single([0.9, 0.8, 0.7, 0.6, 0.0, 0.0])
    obs = adata.obs
    assert [obs[f"cas_cell_type_name_{n}"].iloc[0] for n in range(1, 5)] == [
        "CL:0000004",
        "CL:0000003",
        "CL:0000002",
        "CL:0000001",
    ]
    assert [obs[f"cas_cell_type_label_{n}"].iloc[0] for n in range(1, 5)] == [
        "label_CL:0000004",
        "label_CL:0000003",
        "label_CL:0000002",
        "label_CL:0000001",
    ]
    assert np.allclose([obs[f"cas_cell_type_score_{n}"].iloc[0] for n in range(1, 5)], [0.6, 0.7, 0.8, 0.9])
    assert not [c for c in obs.columns if c.startswith("cas_cell_type_") and c.endswith(("_5", "_0"))]


def test_identical_l3_l4_term_fills_two_positions():
    adata = _run_single([0.0, 0.0, 0.0, 0.0, 0.7, 0.0])
    assert adata.obs["cas_cell_type_name_1"].iloc[0] == "CL:0000005"
    assert adata.obs["cas_cell_type_name_2"].iloc[0] == "CL:0000005"
    assert adata.obs["cas_cell_type_name_3"].iloc[0] == "CL:0000002"


def test_fill_missing_levels_flag_controls_undefined_levels():
    # CL:0000003 selected (depth 3): lookup l4 is undefined.
    row = [0.9, 0.8, 0.7, 0.0, 0.0, 0.0]
    filled = _run_single(row, fill_missing_levels=True).obs
    assert filled["cas_cell_type_name_1"].iloc[0] == "CL:0000003"
    assert filled["cas_cell_type_label_1"].iloc[0] == "label_CL:0000003"
    assert filled["cas_cell_type_score_1"].iloc[0] == 0.7

    unfilled = _run_single(row, fill_missing_levels=False).obs
    assert pd.isna(unfilled["cas_cell_type_name_1"].iloc[0])
    assert pd.isna(unfilled["cas_cell_type_label_1"].iloc[0])
    assert pd.isna(unfilled["cas_cell_type_score_1"].iloc[0])
    assert unfilled["cas_cell_type_name_2"].iloc[0] == "CL:0000003"


def test_no_eligible_candidates_do_not_fabricate_confidence():
    adata = _run_single([0.1] * 6)
    for position in range(1, 5):
        for suffix in ("name", "label", "score"):
            assert pd.isna(adata.obs[f"cas_cell_type_{suffix}_{position}"].iloc[0])


def test_cluster_outputs_use_aggregated_scores():
    cache = CellGranularityCache(TEST_GRANULARITY_RESOURCE)
    cl = _make_cl()

    # Cluster A has mean depth-2 score 0.5; cluster B has mean 0.3.
    scores = np.zeros((4, len(TEST_CL_NAMES)))
    scores[:, 1] = [0.4, 0.6, 0.2, 0.4]
    scores[:, 0] = [0.8, 0.8, 0.8, 0.8]
    adata = _make_adata(scores)
    adata.obs["cluster_label"] = pd.Categorical(["A", "A", "B", "B"])

    compute_most_granular_top_k_calls_cluster(
        adata=adata,
        cl=cl,
        granularity_cache=cache,
        min_acceptable_score=0.2,
        cluster_label_obs_column="cluster_label",
        obs_prefix="cas_cluster_type",
    )

    # CL:0000002 selected; lookup l3/l4 undefined so positions 1-3 repeat it; position 4 is l1.
    for position in (1, 2, 3):
        assert adata.obs[f"cas_cluster_type_name_{position}"].tolist() == ["CL:0000002"] * 4
        assert np.allclose(adata.obs[f"cas_cluster_type_score_{position}"].to_numpy(), [0.5, 0.5, 0.3, 0.3])
    assert adata.obs["cas_cluster_type_name_4"].tolist() == ["CL:0000001"] * 4
    assert np.allclose(adata.obs["cas_cluster_type_score_4"].to_numpy(), [0.8] * 4)


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
        obs_prefix="cas_knn_cell_type",
    )

    expected_score = np.mean(scores[:, 2])
    assert adata.obs["cas_knn_cell_type_name_1"].tolist() == ["CL:0000003"] * 3
    assert np.allclose(pd.to_numeric(adata.obs["cas_knn_cell_type_score_1"]).to_numpy(), expected_score)
    assert adata.obs["cas_knn_cell_type_name_2"].tolist() == ["CL:0000003"] * 3
