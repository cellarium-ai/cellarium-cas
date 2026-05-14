#!/usr/bin/env python3
"""
Build a CAS-compatible ``ontology_response.json`` from Azimuth cell type annotations.

Each cell's Azimuth annotations (one per granularity level) are mapped to CL ontology
terms via the HRA crosswalk and then propagated up the ontology ancestor chain.
Propagation follows these rules:

- Levels are processed **most granular first**.
- For each level, the mapped CL term and all its ancestors receive
  ``score = max(current_score, level_score)``.
- When a less-granular level's CL term is not a confirmed ancestor of the more-granular
  level's CL term (crosswalk linearity violation), a warning is emitted once per unique
  pair and propagation continues independently for both.

The resulting ``CellTypeOntologyAwareResults`` mimics the CAS API response format so
that ``cellarium-cas benchmark ontology-aware`` can evaluate Azimuth on equal footing.
Synthetic fields not produced by Azimuth are set as: ``total_neighbors=1``,
``total_neighbors_unrecognized=0``, ``total_weight=<most_granular_level_score>``.

Usage
-----
::

    from cellarium.cas.benchmarking.azimuth.helpers.build_ontology_response import build_ontology_response

    build_ontology_response(
        inferred_labels_path="/data/azimuth_annotate_dir/inferred_labels.csv",
        ontology_resource_path="/path/to/ontology_resource.json",
        azimuth_ref_name="pbmcref",
        output_dir="/data/azimuth_annotate_dir",
    )

Run :func:`map_azimuth_to_cas_labels` first to produce ``inferred_labels.csv``.
The ontology resource is the ``ontology_resource.json`` saved by
``cellarium-cas annotate --save-ontology-resource``.
"""

import json
import typing as t
import warnings
from pathlib import Path

from cellarium.cas.models import CellTypeOntologyAwareResults
from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache


def _build_parents_map(children_dict: t.Dict[str, t.List[str]]) -> t.Dict[str, t.List[str]]:
    """Invert ``children_dict`` (child → parents) for ancestor traversal."""
    parents_map: t.Dict[str, t.List[str]] = {}
    for parent, children in children_dict.items():
        for child in children:
            parents_map.setdefault(child, []).append(parent)
    return parents_map


def _get_node_and_ancestors(cl_id: str, parents_map: t.Dict[str, t.List[str]]) -> t.Set[str]:
    """Return *cl_id* and all its ontology ancestors via BFS."""
    visited: t.Set[str] = {cl_id}
    queue = [cl_id]
    while queue:
        node = queue.pop(0)
        for parent in parents_map.get(node, []):
            if parent not in visited:
                visited.add(parent)
                queue.append(parent)
    return visited


def build_ontology_response(
    inferred_labels_path: str,
    ontology_resource_path: str,
    azimuth_ref_name: str,
    output_dir: str,
) -> CellTypeOntologyAwareResults:
    """
    Build a CAS-compatible ``CellTypeOntologyAwareResults`` from an ``inferred_labels.csv``
    produced by :func:`map_azimuth_to_cas_labels`.

    Reads ``cas_cell_type_label_k`` and ``cas_cell_type_score_k`` columns (rank 1 = most
    granular) and propagates each mapped CL term's score to all its ontology ancestors
    using ``max``-merge, so ancestor nodes reflect the highest-confidence prediction that
    flows through them.

    :param inferred_labels_path: Path to ``inferred_labels.csv`` from
        :func:`map_azimuth_to_cas_labels` (row index = barcodes, columns
        ``cas_cell_type_label_k`` / ``cas_cell_type_score_k``).
    :param ontology_resource_path: Path to ``ontology_resource.json`` (saved by
        ``cellarium-cas annotate --save-ontology-resource``).  Loaded to build the
        :class:`~cellarium.cas.postprocessing.cell_ontology.CellOntologyCache` and
        copied as-is to *output_dir*.
    :param azimuth_ref_name: Azimuth reference name (e.g. ``"pbmcref"``), written into
        ``model_name`` as ``azimuth_<ref_name>``.
    :param output_dir: Directory to write ``ontology_response.json`` and
        ``ontology_resource.json`` into (created if absent).

    :returns: The constructed :class:`~cellarium.cas.models.CellTypeOntologyAwareResults`.
    :raises ValueError: If no ``cas_cell_type_label_k`` columns are found.
    """
    import re

    import pandas as pd

    # --- load ontology resource and build cache ---
    with open(ontology_resource_path) as f:
        ontology_resource = json.load(f)
    cell_ontology_cache = CellOntologyCache(ontology_resource)

    # --- load inferred labels ---
    labels_df = pd.read_csv(inferred_labels_path, index_col=0)
    labels_df.index = labels_df.index.astype(str)

    top_k = max(
        (int(m.group(1)) for col in labels_df.columns if (m := re.match(r"^cas_cell_type_label_(\d+)$", col))),
        default=0,
    )
    if top_k == 0:
        raise ValueError(
            f"No 'cas_cell_type_label_k' columns found in {inferred_labels_path}. "
            "Run map_azimuth_to_cas_labels.py first."
        )

    # --- pre-build ancestor traversal structures ---
    parents_map = _build_parents_map(cell_ontology_cache.children_dict)
    known_cl_ids = set(cell_ontology_cache.cl_names_to_idx_map.keys())
    ancestors_cache: t.Dict[str, t.Set[str]] = {}

    def get_ancestors_cached(cl_id: str) -> t.Set[str]:
        if cl_id not in ancestors_cache:
            ancestors_cache[cl_id] = _get_node_and_ancestors(cl_id, parents_map)
        return ancestors_cache[cl_id]

    # --- per-cell processing ---
    warned_missing_cl_ids: t.Set[str] = set()
    warned_nonlinear_pairs: t.Set[t.Tuple[str, str]] = set()

    annotations: t.List[CellTypeOntologyAwareResults.OntologyAwareAnnotation] = []

    for barcode, row in labels_df.iterrows():
        # Collect (cl_id, score) per rank; rank 1 = most granular
        level_cl_ids: t.List[t.Tuple[t.Optional[str], float]] = []
        for k in range(1, top_k + 1):
            raw_cl_id = row.get(f"cas_cell_type_label_{k}")
            raw_score = row.get(f"cas_cell_type_score_{k}")

            score = float(raw_score) if raw_score is not None and not pd.isna(raw_score) else 0.0
            cl_id: t.Optional[str] = None

            if raw_cl_id is not None and not pd.isna(raw_cl_id):
                cl_id = str(raw_cl_id)
                if cl_id not in known_cl_ids:
                    if cl_id not in warned_missing_cl_ids:
                        warnings.warn(
                            f"CL ID '{cl_id}' (rank {k}) is not present in the cell ontology cache — skipped."
                        )
                        warned_missing_cl_ids.add(cl_id)
                    cl_id = None

            level_cl_ids.append((cl_id, score))

        # Linearity check (warn once per unique pair, A+C strategy)
        valid_ids = [(cl_id, score) for cl_id, score in level_cl_ids if cl_id is not None]
        for i in range(len(valid_ids) - 1):
            more_granular_cl = valid_ids[i][0]
            less_granular_cl = valid_ids[i + 1][0]
            pair = (more_granular_cl, less_granular_cl)
            if pair not in warned_nonlinear_pairs:
                ancestors_of_more_granular = get_ancestors_cached(more_granular_cl)
                if less_granular_cl not in ancestors_of_more_granular:
                    warnings.warn(
                        f"Crosswalk hierarchy violation: '{less_granular_cl}' is not an ancestor of "
                        f"'{more_granular_cl}' in the cell ontology — propagating independently."
                    )
                    warned_nonlinear_pairs.add(pair)

        # Score propagation: most-granular first, max-merge at shared ancestors
        node_scores: t.Dict[str, float] = {}
        for cl_id, score in level_cl_ids:
            if cl_id is None:
                continue
            for ancestor in get_ancestors_cached(cl_id):
                if ancestor in known_cl_ids:
                    node_scores[ancestor] = max(node_scores.get(ancestor, 0.0), score)

        matches = [
            CellTypeOntologyAwareResults.Match(
                cell_type_ontology_term_id=cl_id,
                cell_type=cell_ontology_cache.cl_names_to_labels_map.get(cl_id, cl_id),
                score=score,
            )
            for cl_id, score in node_scores.items()
        ]

        # total_weight: score of the most-granular non-null level
        most_granular_score = next((score for cl_id, score in level_cl_ids if cl_id is not None), 0.0)

        annotations.append(
            CellTypeOntologyAwareResults.OntologyAwareAnnotation(
                query_cell_id=barcode,
                matches=matches,
                total_weight=most_granular_score,
                total_neighbors=1,
                total_neighbors_unrecognized=0,
            )
        )

    result = CellTypeOntologyAwareResults(
        data=annotations,
        model_name=f"azimuth_{azimuth_ref_name}",
    )

    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    with open(output_dir_path / "ontology_response.json", "w") as f:
        f.write(result.model_dump_json())
    with open(output_dir_path / "ontology_resource.json", "w") as f:
        json.dump(ontology_resource, f)

    return result
