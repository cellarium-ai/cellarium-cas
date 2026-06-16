"""
Shared test hierarchy for confusion-matrix and F-measure unit tests.

Hierarchy:
    CL:0000000 (root -- excluded from LABELS and all ancestor sets)
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

from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache

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
