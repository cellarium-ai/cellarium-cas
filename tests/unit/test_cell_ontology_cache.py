"""
Unit tests for CellOntologyCache and its internal BFS helper.

Hierarchy used throughout:

    CL:0000000 (root)
        └── CL:0000001 (parent)
                ├── CL:0000002 (child_a)  -- leaf
                └── CL:0000003 (child_b)
                        └── CL:0000004 (grandchild)  -- leaf

Ancestor sets (self-inclusive, root-inclusive):
    CL:0000000 -> {CL:0000000}
    CL:0000001 -> {CL:0000000, CL:0000001}
    CL:0000002 -> {CL:0000000, CL:0000001, CL:0000002}
    CL:0000003 -> {CL:0000000, CL:0000001, CL:0000003}
    CL:0000004 -> {CL:0000000, CL:0000001, CL:0000003, CL:0000004}

Descendant sets (self-inclusive):
    CL:0000000 -> {CL:0000000, CL:0000001, CL:0000002, CL:0000003, CL:0000004}
    CL:0000001 -> {CL:0000001, CL:0000002, CL:0000003, CL:0000004}
    CL:0000002 -> {CL:0000002}
    CL:0000003 -> {CL:0000003, CL:0000004}
    CL:0000004 -> {CL:0000004}
"""

import typing as t

from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache, _bfs_reachable

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


# ---------------------------------------------------------------------------
# _bfs_reachable (unit)
# ---------------------------------------------------------------------------


def test_bfs_isolated_node_returns_singleton():
    assert _bfs_reachable("A", {}) == frozenset({"A"})


def test_bfs_node_absent_from_adjacency_returns_singleton():
    assert _bfs_reachable("A", {"B": ["C"]}) == frozenset({"A"})


def test_bfs_linear_chain():
    adj = {"A": ["B"], "B": ["C"]}
    assert _bfs_reachable("A", adj) == frozenset({"A", "B", "C"})


def test_bfs_start_from_middle_does_not_backtrack():
    adj = {"A": ["B"], "B": ["C"]}
    assert _bfs_reachable("B", adj) == frozenset({"B", "C"})


def test_bfs_diamond_dag_visits_shared_node_once():
    # A -> B, A -> C, B -> D, C -> D
    adj = {"A": ["B", "C"], "B": ["D"], "C": ["D"]}
    assert _bfs_reachable("A", adj) == frozenset({"A", "B", "C", "D"})


def test_bfs_cycle_does_not_loop_forever():
    adj = {"A": ["B"], "B": ["A"]}
    assert _bfs_reachable("A", adj) == frozenset({"A", "B"})


def test_bfs_result_is_frozenset():
    assert isinstance(_bfs_reachable("X", {}), frozenset)


# ---------------------------------------------------------------------------
# CellOntologyCache — ancestor building
# ---------------------------------------------------------------------------


def test_ancestors_raw_includes_root():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    for term in MOCK_ONTOLOGY_RESOURCE["cl_names"]:
        if term == "CL:0000000":
            continue
        assert "CL:0000000" in cache.ancestors_dict[term], f"Root missing from raw ancestors of {term}"


def test_ancestors_raw_includes_self():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    for term in MOCK_ONTOLOGY_RESOURCE["cl_names"]:
        assert term in cache.ancestors_dict[term], f"{term} not in its own ancestor set"


def test_ancestors_remove_root_flag():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    for term in MOCK_ONTOLOGY_RESOURCE["cl_names"]:
        assert "CL:0000000" not in cache.get_ancestors(term, remove_root=True)


def test_ancestors_remove_self_flag():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    for term in MOCK_ONTOLOGY_RESOURCE["cl_names"]:
        assert term not in cache.get_ancestors(term, remove_self=True)


def test_ancestors_grandchild_remove_root():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    assert cache.get_ancestors("CL:0000004", remove_root=True) == frozenset({"CL:0000001", "CL:0000003", "CL:0000004"})


def test_ancestors_child_a_remove_root():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    assert cache.get_ancestors("CL:0000002", remove_root=True) == frozenset({"CL:0000001", "CL:0000002"})


def test_ancestors_unknown_term_fallback():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    assert cache.get_ancestors("CL:9999999") == frozenset({"CL:9999999"})


# ---------------------------------------------------------------------------
# CellOntologyCache — descendant building
# ---------------------------------------------------------------------------


def test_descendants_raw_includes_self():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    for term in MOCK_ONTOLOGY_RESOURCE["cl_names"]:
        assert term in cache.descendants_dict[term], f"{term} not in its own descendant set"


def test_descendants_leaf_nodes_are_singletons():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    for leaf in ("CL:0000002", "CL:0000004"):
        assert cache.descendants_dict[leaf] == frozenset({leaf})


def test_descendants_root_contains_all_terms():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    assert cache.descendants_dict["CL:0000000"] == frozenset(MOCK_ONTOLOGY_RESOURCE["cl_names"])


def test_descendants_parent_node():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    assert cache.descendants_dict["CL:0000001"] == frozenset({"CL:0000001", "CL:0000002", "CL:0000003", "CL:0000004"})


def test_descendants_child_b():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    assert cache.get_descendants("CL:0000003") == frozenset({"CL:0000003", "CL:0000004"})


def test_descendants_remove_self_flag():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    for term in MOCK_ONTOLOGY_RESOURCE["cl_names"]:
        assert term not in cache.get_descendants(term, remove_self=True)


def test_descendants_remove_root_flag():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    for term in MOCK_ONTOLOGY_RESOURCE["cl_names"]:
        assert "CL:0000000" not in cache.get_descendants(term, remove_root=True)


def test_descendants_unknown_term_fallback():
    cache = CellOntologyCache(MOCK_ONTOLOGY_RESOURCE)
    assert cache.get_descendants("CL:9999999") == frozenset({"CL:9999999"})
