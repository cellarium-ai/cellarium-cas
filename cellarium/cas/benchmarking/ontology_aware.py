import typing as t

import pandas as pd

from cellarium.cas.models import CellTypeOntologyAwareResults


class _OntologyGraph(t.NamedTuple):
    """Precomputed ontology graph structures used for hop-neighbourhood evaluation."""

    children: t.Dict[str, t.Set[str]]
    parents: t.Dict[str, t.Set[str]]
    term_ancestors: t.Dict[str, t.FrozenSet[str]]
    term_descendants: t.Dict[str, t.FrozenSet[str]]
    all_terms: t.FrozenSet[str]


def _reachable(start: str, adjacency: t.Dict[str, t.Set[str]]) -> t.FrozenSet[str]:
    """BFS from *start* through *adjacency*. Includes *start* (a node is its own ancestor/descendant)."""
    visited: t.Set[str] = {start}
    frontier: t.Set[str] = {start}
    while frontier:
        next_frontier: t.Set[str] = set()
        for node in frontier:
            for neighbor in adjacency.get(node, ()):
                if neighbor not in visited:
                    visited.add(neighbor)
                    next_frontier.add(neighbor)
        frontier = next_frontier
    return frozenset(visited)


def _precompute_graph(ontology_resource: t.Dict[str, t.Any]) -> _OntologyGraph:
    """
    Build parent/child adjacency and precompute ancestor/descendant sets for all ontology terms.

    :param ontology_resource: Raw cell ontology resource dict with keys ``cl_names`` and
        ``children_dictionary`` (same shape as ``CellOntologyCache.__init__`` accepts).

    :return: :class:`_OntologyGraph` namedtuple with precomputed graph structures.
    """
    children: t.Dict[str, t.Set[str]] = {k: set(v) for k, v in ontology_resource["children_dictionary"].items()}

    parents: t.Dict[str, t.Set[str]] = {}
    for parent_term, child_terms in children.items():
        for child in child_terms:
            parents.setdefault(child, set()).add(parent_term)

    all_terms: t.FrozenSet[str] = frozenset(ontology_resource["cl_names"])

    term_ancestors: t.Dict[str, t.FrozenSet[str]] = {term: _reachable(term, parents) for term in all_terms}
    term_descendants: t.Dict[str, t.FrozenSet[str]] = {term: _reachable(term, children) for term in all_terms}

    return _OntologyGraph(
        children=children,
        parents=parents,
        term_ancestors=term_ancestors,
        term_descendants=term_descendants,
        all_terms=all_terms,
    )


def _build_hop_neighborhoods(
    gt_term: str,
    num_hops: int,
    graph: _OntologyGraph,
) -> t.List[t.Dict[str, t.FrozenSet[str]]]:
    """
    Build cumulative bidirectional hop-N neighbourhoods around *gt_term* for N in 0..num_hops.

    Hop-0 contains only *gt_term* itself (exact match). Hop-N expands by one bidirectional
    step from the previous frontier.  Each entry is a dict with:

    - ``nodes``: all ontology terms reachable within N bidirectional steps of *gt_term*
    - ``all_ancestors``: union of ancestor sets (self-inclusive) of all nodes in this hop
    - ``all_descendants``: union of descendant sets (self-inclusive) of all nodes in this hop

    :param gt_term: Ground truth CL term ID.
    :param num_hops: Maximum hop distance to build (inclusive).
    :param graph: Precomputed :class:`_OntologyGraph`.

    :return: List of length ``num_hops + 1`` (index N = hop N).
    """
    hops = []
    visited: t.Set[str] = {gt_term}
    frontier: t.Set[str] = {gt_term}

    for n in range(num_hops + 1):
        if n > 0:
            next_frontier: t.Set[str] = set()
            for node in frontier:
                next_frontier.update(graph.children.get(node, ()))
                next_frontier.update(graph.parents.get(node, ()))
            frontier = next_frontier - visited
            visited |= frontier

        nodes = frozenset(visited)
        all_ancestors = frozenset(ancestor for node in nodes for ancestor in graph.term_ancestors.get(node, ()))
        all_descendants = frozenset(descendant for node in nodes for descendant in graph.term_descendants.get(node, ()))
        hops.append({"nodes": nodes, "all_ancestors": all_ancestors, "all_descendants": all_descendants})

    return hops


def _precision(tp: float, fp: float) -> float:
    return tp / (tp + fp) if (tp + fp) > 0.0 else 0.0


def _f1(precision: float, recall: float) -> float:
    return (2.0 * precision * recall) / (precision + recall) if (precision + recall) > 0.0 else 0.0


def _compute_cell_metrics(
    matches: t.List[t.Any],
    hops: t.List[t.Dict[str, t.FrozenSet[str]]],
    graph: _OntologyGraph,
    num_hops: int,
) -> t.Dict[str, float]:
    """
    Compute per-hop sensitivity, specificity, and F1 for a single cell's prediction matches.

    At hop-N a predicted term is:

    - **True positive**: the predicted term falls within the hop-N neighbourhood.
    - **False positive**: the predicted term lies entirely outside the neighbourhood's
      ancestor/descendant lineage *and* its own descendant subtree shares no overlap with
      the neighbourhood's descendants.
    - **Neither** (related but outside neighbourhood): skipped for TP/FP counting.

    When multiple matches exist, the maximum score among TPs (or FPs) is used.

    :param matches: List of :class:`CellTypeOntologyAwareResults.Match` objects for one cell.
    :param hops: Precomputed hop neighbourhoods from :func:`_build_hop_neighborhoods`.
    :param graph: Precomputed :class:`_OntologyGraph`.
    :param num_hops: Maximum hop index to evaluate.

    :return: Dict with keys ``hop_{N}_sensitivity``, ``hop_{N}_specificity``, ``hop_{N}_f1_score``
        for N in 0..num_hops.
    """
    true_positives = [0.0] * (num_hops + 1)
    false_positives = [0.0] * (num_hops + 1)

    for match in matches:
        match_term = match.cell_type_ontology_term_id
        match_score = match.score

        if match_term not in graph.term_ancestors:
            continue

        match_all_descendants = graph.term_descendants[match_term]

        for i, hop in enumerate(hops):
            if match_term in hop["nodes"]:
                true_positives[i] = max(match_score, true_positives[i])
            elif (
                match_term not in hop["all_ancestors"]
                and match_term not in hop["all_descendants"]
                and not match_all_descendants.intersection(hop["all_descendants"])
            ):
                false_positives[i] = max(match_score, false_positives[i])

    result: t.Dict[str, float] = {}
    for i in range(num_hops + 1):
        tp = true_positives[i]
        fp = false_positives[i]
        prec = _precision(tp, fp)
        result[f"hop_{i}_sensitivity"] = tp
        result[f"hop_{i}_specificity"] = 1.0 - fp
        result[f"hop_{i}_f1_score"] = _f1(prec, tp)

    return result


def compute_ontology_aware_metrics(
    response: CellTypeOntologyAwareResults,
    ground_truths: t.List[str],
    ontology_resource: t.Dict[str, t.Any],
    num_hops: int = 4,
    cell_level: bool = False,
) -> pd.DataFrame:
    """
    Compute ontology-aware benchmarking metrics from a CAS ontology-aware response.

    For each cell, evaluates predictions at hop levels 0 through ``num_hops``. At hop N,
    a prediction is a true positive if the predicted cell type falls within N bidirectional
    steps of the ground truth in the Cell Ontology graph. A prediction is a false positive
    only if it lies entirely outside the hop neighbourhood's ancestor/descendant lineage
    *and* its own descendant subtree shares no overlap with the neighbourhood's descendants.

    :param response: CAS ontology-aware response object (``CellTypeOntologyAwareResults``).
    :param ground_truths: Ground truth CL term IDs (e.g. ``"CL:0000121"``), positionally aligned
        with ``response.data``. Must have the same length as ``response.data``.
    :param ontology_resource: Raw cell ontology resource dict with at minimum keys
        ``cl_names`` (list of CL term IDs) and ``children_dictionary`` (term ID to list of child IDs).
        This is the same dict accepted by ``CellOntologyCache.__init__``.
    :param num_hops: Maximum hop distance to evaluate (default 4).
    :param cell_level: If ``False`` (default), return a summary DataFrame with mean sensitivity,
        specificity, and F1 per hop. If ``True``, return a per-cell DataFrame with
        per-hop metrics for each cell.

    :return:
        - Summary (``cell_level=False``): DataFrame with columns
          ``hop_{N}_sensitivity``, ``hop_{N}_specificity``, ``hop_{N}_f1_score`` for N in 0..num_hops,
          plus ``n_cells``.
        - Cell-level (``cell_level=True``): DataFrame with columns ``query_cell_id``,
          ``ground_truth``, and per-hop metric columns.

    :raises ValueError: If ``len(ground_truths) != len(response.data)`` or if any ground truth
        term is not present in the ontology resource.
    """
    if len(ground_truths) != len(response.data):
        raise ValueError(
            f"Length mismatch: ground_truths has {len(ground_truths)} entries "
            f"but response.data has {len(response.data)} entries."
        )

    all_terms = frozenset(ontology_resource["cl_names"])
    unrecognized = [gt for gt in ground_truths if gt not in all_terms]
    if unrecognized:
        raise ValueError(
            f"The following ground truth terms are not present in the ontology resource: {sorted(set(unrecognized))}. "
            f"Ensure all ground truth labels are valid CL term IDs."
        )

    graph = _precompute_graph(ontology_resource)
    hop_cache: t.Dict[str, t.List[t.Dict[str, t.FrozenSet[str]]]] = {}
    cell_rows: t.List[t.Dict[str, t.Any]] = []

    for annotation, gt_term in zip(response.data, ground_truths):
        if gt_term not in hop_cache:
            hop_cache[gt_term] = _build_hop_neighborhoods(gt_term, num_hops, graph)

        hops = hop_cache[gt_term]
        cell_metrics = _compute_cell_metrics(annotation.matches, hops, graph, num_hops)

        row: t.Dict[str, t.Any] = {"query_cell_id": annotation.query_cell_id, "ground_truth": gt_term}
        row.update(cell_metrics)
        cell_rows.append(row)

    if not cell_rows:
        return pd.DataFrame()

    cell_df = pd.DataFrame(cell_rows)

    if cell_level:
        return cell_df

    metric_cols = [
        f"hop_{i}_{metric}" for i in range(num_hops + 1) for metric in ("sensitivity", "specificity", "f1_score")
    ]
    summary = cell_df[metric_cols].mean().to_frame().T
    summary.insert(0, "n_cells", len(cell_df))
    summary = summary.reset_index(drop=True)
    return summary
