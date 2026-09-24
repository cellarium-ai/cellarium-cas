import typing as t
from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from enum import Enum
from operator import itemgetter

import numpy as np
import scipy.sparse as sp
from anndata import AnnData

from cellarium.cas.logging import logger
from cellarium.cas.models import CellTypeOntologyAwareResults

from .cell_ontology.cell_ontology_cache import CL_CELL_ROOT_NODE, CellOntologyCache
from .common import get_obs_indices_for_cluster

# AnnData-related constants
CAS_CL_SCORES_ANNDATA_OBSM_KEY = "cas_cl_scores"
CAS_METADATA_ANNDATA_UNS_KEY = "cas_metadata"

# PhyloXML-related constants
CAS_SCORE_PROPERTY_REF = "CAS:score"
CAS_FRACTION_PROPERTY_REF = "CAS:fraction"
CAS_CL_LABEL_PROPERTY_REF = "CAS:cl_label"


def convert_cas_ontology_aware_response_to_score_matrix(
    adata: AnnData, cas_ontology_aware_response: CellTypeOntologyAwareResults, cl: CellOntologyCache
) -> sp.csr_matrix:
    """
    Generate a sparse matrix of CAS ontology-aware scores.

    This function takes an AnnData object, a list of CAS ontology-aware responses, and a CellOntologyCache object
    and generates a sparse matrix of CAS ontology-aware scores. The sparse matrix represents the evidence scores
    of all ontology terms (columns) for each cell (rows).

    :param adata: An AnnData object containing the query cells.
    :type adata: AnnData

    :param cas_ontology_aware_response: A list of CAS ontology-aware responses.
    :type cas_ontology_aware_response: CellTypeOntologyAwareResults

    :param cl: A CellOntologyCache object containing the cell ontology information.
    :type cl: CellOntologyCache

    :return: A sparse matrix of CAS ontology-aware scores.
    :rtype: sp.csr_matrix
    """
    row = []
    col = []
    data = []

    obs_values = adata.obs.index.values
    for obs_idx, cas_cell_response in enumerate(cas_ontology_aware_response.data):
        assert cas_cell_response.query_cell_id == obs_values[obs_idx]
        for match in cas_cell_response.matches:
            row.append(obs_idx)
            col.append(cl.cl_names_to_idx_map[match.cell_type_ontology_term_id])
            data.append(match.score)

    n_obs = len(cas_ontology_aware_response.data)
    n_cl_names = len(cl.cl_names)
    return sp.coo_matrix((data, (row, col)), shape=(n_obs, n_cl_names)).tocsr()


def insert_cas_ontology_aware_response_into_adata(
    cas_ontology_aware_response: CellTypeOntologyAwareResults,
    adata: AnnData,
    cl: CellOntologyCache,
    ontology_resource_name: t.Optional[str] = None,
) -> None:
    """
    Inserts Cellarium CAS ontology aware response into `obsm` property of a provided AnnData file as a
    :class:`scipy.sparse.csr_matrix` named `cas_cl_scores`.

    :param cas_ontology_aware_response: The Cellarium CAS ontology aware response.
    :type cas_ontology_aware_response: CellTypeOntologyAwareResults

    :param adata: The AnnData object to insert the response into.
    :type adata: AnnData

    :param cl: The CellOntologyCache object containing cell ontology term names and labels.
    :type cl: CellOntologyCache

    :param ontology_resource_name: Optional name of the ontology resource, stored in adata.uns for later use.
    :type ontology_resource_name: str, optional

    :return: None
    """
    adata.obsm[CAS_CL_SCORES_ANNDATA_OBSM_KEY] = convert_cas_ontology_aware_response_to_score_matrix(
        adata, cas_ontology_aware_response, cl
    )
    cas_metadata = {"cl_names": cl.cl_names, "cl_labels": cl.cl_labels}
    if ontology_resource_name is not None:
        cas_metadata["ontology_resource_name"] = ontology_resource_name
    adata.uns[CAS_METADATA_ANNDATA_UNS_KEY] = cas_metadata


class CellOntologyScoresAggregationOp(Enum):
    """The aggregation operation to apply to the scores of each ontology term."""

    MEAN = 1
    MEDIAN = 2


class CellOntologyScoresAggregationDomain(Enum):
    """The domain over which to aggregate the scores of each ontology term."""

    ALL_CELLS = 1
    OVER_THRESHOLD = 2


@dataclass
class AggregatedCellOntologyScores:
    """The aggregated ontology scores and the fraction of cells that have non-zero scores for each term."""

    aggregation_op: CellOntologyScoresAggregationOp
    aggregation_domain: CellOntologyScoresAggregationDomain
    threshold: float
    aggregated_scores_c: np.ndarray
    fraction_over_threshold_c: np.ndarray
    cl_labels: list
    cl_names: list


def get_aggregated_cas_ontology_aware_scores(
    adata: AnnData,
    obs_indices: t.Optional[t.Sequence],
    aggregation_op: CellOntologyScoresAggregationOp = CellOntologyScoresAggregationOp.MEAN,
    aggregation_domain: CellOntologyScoresAggregationDomain = CellOntologyScoresAggregationDomain.ALL_CELLS,
    threshold: float = 1e-4,
) -> AggregatedCellOntologyScores:
    """
    Aggregates Cellarium CAS ontology-aware scores over a provided group of cells.

    :param adata: Annotated data object containing CAS scores.
    :param obs_indices: Indices of the cells to consider for aggregation. If None, all cells are considered.
    :param aggregation_op: The aggregation operation to apply to the CAS scores.
      Default is CellOntologyScoresAggregationOp.MEAN.
    :param aggregation_domain: The domain over which to perform the aggregation.
      Default is CellOntologyScoresAggregationDomain.ALL_CELLS.
    :param threshold: The threshold value for considering a CAS score as non-zero.
      Default is 1e-4.
    :return: An instance of AggregatedCellOntologyScores.
    """
    assert CAS_CL_SCORES_ANNDATA_OBSM_KEY in adata.obsm

    if obs_indices is not None:
        assert len(obs_indices) > 0
        sliced_cas_scores_nc = adata.obsm[CAS_CL_SCORES_ANNDATA_OBSM_KEY][obs_indices]
        n_selected_cells = len(obs_indices)
    else:  # all cells
        sliced_cas_scores_nc = adata.obsm[CAS_CL_SCORES_ANNDATA_OBSM_KEY]
        n_selected_cells = adata.n_obs

    mask_aggregation_nc = (sliced_cas_scores_nc > threshold).toarray()
    fraction_over_threshold_c = mask_aggregation_nc.sum(0) / n_selected_cells

    if aggregation_op == CellOntologyScoresAggregationOp.MEAN:
        aggregation_op_func = np.mean
    elif aggregation_op == CellOntologyScoresAggregationOp.MEDIAN:
        aggregation_op_func = np.median
    else:
        raise ValueError

    sliced_cas_scores_dense_nc = sliced_cas_scores_nc.toarray()
    n_cols = sliced_cas_scores_dense_nc.shape[1]
    aggregated_scores_c = np.zeros((n_cols,))
    for c in range(n_cols):
        if aggregation_domain == CellOntologyScoresAggregationDomain.ALL_CELLS:
            value = aggregation_op_func(sliced_cas_scores_dense_nc[:, c])
        elif aggregation_domain == CellOntologyScoresAggregationDomain.OVER_THRESHOLD:
            data = sliced_cas_scores_dense_nc[:, c][mask_aggregation_nc[:, c]]
            value = aggregation_op_func(data) if data.size else 0.0
        else:
            raise ValueError
        aggregated_scores_c[c] = value

    return AggregatedCellOntologyScores(
        aggregation_op=aggregation_op,
        aggregation_domain=aggregation_domain,
        threshold=threshold,
        aggregated_scores_c=aggregated_scores_c,
        fraction_over_threshold_c=fraction_over_threshold_c,
        cl_labels=adata.uns[CAS_METADATA_ANNDATA_UNS_KEY]["cl_labels"],
        cl_names=adata.uns[CAS_METADATA_ANNDATA_UNS_KEY]["cl_names"],
    )


def convert_aggregated_cell_ontology_scores_to_rooted_tree(
    aggregated_scores: AggregatedCellOntologyScores,
    cl: CellOntologyCache,
    root_cl_name: str = CL_CELL_ROOT_NODE,
    min_fraction: float = 0.0,
    hidden_cl_names_set: t.Optional[t.Set[str]] = None,
) -> OrderedDict:
    """
    Convert aggregated cell ontology scores to a rooted tree.

    :param aggregated_scores: The aggregated cell ontology scores.
    :type aggregated_scores: AggregatedCellOntologyScores
    :param cl: The cell ontology cache.
    :type cl: CellOntologyCache
    :param root_cl_name: The name of the root cell ontology.
    :type root_cl_name: str
    :param min_fraction: The minimum fraction threshold, defaults to 0.0.
    :type min_fraction: float, optional
    :param hidden_cl_names_set: A set of hidden cell ontology names, defaults to an empty set.
    :type hidden_cl_names_set: set[str], optional
    :return: The rooted tree structure as an ordered dictionary.
    :rtype: OrderedDict
    """

    if hidden_cl_names_set is None:
        hidden_cl_names_set = set()

    tree_dict = OrderedDict()
    score_dict = {
        cl_name: score for cl_name, score in zip(aggregated_scores.cl_names, aggregated_scores.aggregated_scores_c)
    }
    fraction_dict = {
        cl_name: score
        for cl_name, score in zip(aggregated_scores.cl_names, aggregated_scores.fraction_over_threshold_c)
    }

    def build_subtree(node_dict: OrderedDict, node_name: str) -> OrderedDict:
        if node_name in hidden_cl_names_set:
            return node_dict
        if score_dict[node_name] <= aggregated_scores.threshold:
            return node_dict
        if fraction_dict[node_name] <= min_fraction:
            return node_dict
        node_dict[node_name] = {"score": score_dict[node_name], "fraction": fraction_dict[node_name]}
        children_nodes = cl.children_dict.get(node_name, [])
        if len(children_nodes) > 0:
            node_dict[node_name]["children"] = OrderedDict()
            for children_node_name in children_nodes:
                build_subtree(node_dict[node_name]["children"], children_node_name)
        return node_dict

    tree_dict = build_subtree(tree_dict, root_cl_name)
    # Validate that this is actually a rooted tree and if not recalculate with the base cell node
    if len(tree_dict) == 1:  # singly-rooted tree
        return tree_dict
    else:
        raise ValueError("The tree is not singly-rooted.")


def generate_phyloxml_from_scored_cell_ontology_tree(
    tree_dict: OrderedDict, tree_name: str, cl: CellOntologyCache, indent: int = 3
) -> str:
    """
    Generate a PhyloXML string representation of a scored cell ontology tree.

    :param tree_dict: The dictionary representation of the scored cell ontology tree.
    :type tree_dict: OrderedDict
    :param tree_name: The name of the tree.
    :type tree_name: str
    :param cl: The CellOntologyCache object containing the cell ontology information.
    :type cl: CellOntologyCache
    :param indent: The number of spaces to use for indentation. Defaults to 3.
    :type indent: int, optional
    :return: The PhyloXML string representation of the scored cell ontology tree.
    :rtype: str
    """

    assert indent > 0
    assert len(tree_dict) == 1  # singly-rooted tree
    root_node_name = list(tree_dict.keys())[0]

    def ind(level: int):
        return " " * indent * level

    phyloxml_header = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.00/phyloxml.xsd" xmlns="http://www.phyloxml.org">\n'
        '<phylogeny rooted="true">\n'
        f"{ind(1)}<name>{tree_name}</name>\n"
    )

    phyloxml_footer = "</phylogeny>\n</phyloxml>\n"

    def _get_subtree_phyloxml_string(subtree_dict: OrderedDict, node_name: str, level: int) -> str:
        phyloxml_string = f"{ind(level)}<clade>\n" f"{ind(level + 1)}<name>{node_name}</name>\n"
        if level > 1:
            phyloxml_string += f"{ind(level + 1)}<branch_length>1.0</branch_length>\n"
        phyloxml_string += (
            f'{ind(level + 1)}<property datatype="xsd:string" ref="{CAS_CL_LABEL_PROPERTY_REF}" applies_to="clade">{cl.cl_names_to_labels_map[node_name]}</property>\n'
            f'{ind(level + 1)}<property datatype="xsd:float" ref="{CAS_SCORE_PROPERTY_REF}" applies_to="clade">{subtree_dict[node_name]["score"]}</property>\n'
            f'{ind(level + 1)}<property datatype="xsd:float" ref="{CAS_FRACTION_PROPERTY_REF}" applies_to="clade">{subtree_dict[node_name]["fraction"]}</property>\n'
        )
        if "children" in subtree_dict[node_name]:
            for child_node_name in subtree_dict[node_name]["children"].keys():
                phyloxml_string += _get_subtree_phyloxml_string(
                    subtree_dict[node_name]["children"], child_node_name, level + 1
                )
        phyloxml_string += f"{ind(level)}</clade>\n"
        return phyloxml_string

    phyloxml_string = phyloxml_header + _get_subtree_phyloxml_string(tree_dict, root_node_name, 1) + phyloxml_footer
    return phyloxml_string


def get_most_granular_top_k_calls(
    aggregated_scores: AggregatedCellOntologyScores,
    cl: CellOntologyCache,
    min_acceptable_score: float,
    top_k: int = 1,
    root_note: str = CL_CELL_ROOT_NODE,
    use_shortest_path: bool = True,
) -> t.List[tuple]:
    if use_shortest_path:
        depth_map = cl.get_shortest_path_lengths_from_target(root_note)
    else:
        depth_map = cl.get_longest_path_lengths_from_target(root_note)
    trunc_list = sorted(
        (
            (score, depth_map[cl_name], cl_name)
            for score, cl_name in zip(aggregated_scores.aggregated_scores_c, aggregated_scores.cl_names)
            if score >= min_acceptable_score
        ),
        key=itemgetter(1, 0),
        reverse=True,
    )[:top_k]
    trunc_list += [(1.0, 0, root_note)] * (top_k - len(trunc_list))
    return trunc_list


def compute_most_granular_top_k_calls_single(
    adata: AnnData,
    cl: CellOntologyCache,
    min_acceptable_score: float,
    top_k: int = 3,
    obs_prefix: str = "cas_cell_type",
    root_note: str = CL_CELL_ROOT_NODE,
    use_shortest_path: bool = True,
):
    top_k_calls_dict = defaultdict(list)
    scores_array_nc = adata.obsm[CAS_CL_SCORES_ANNDATA_OBSM_KEY].toarray()

    # make a wrapped aggregated scores dataclass for the single cell in question
    aggregated_scores = AggregatedCellOntologyScores(
        aggregation_op=CellOntologyScoresAggregationOp.MEAN,
        aggregation_domain=CellOntologyScoresAggregationDomain.ALL_CELLS,
        threshold=0.0,
        aggregated_scores_c=np.zeros((len(cl.cl_names))),
        fraction_over_threshold_c=np.ones((len(cl.cl_names))),
        cl_labels=cl.cl_labels,
        cl_names=cl.cl_names,
    )

    for i_cell in range(adata.n_obs):
        aggregated_scores.aggregated_scores_c = scores_array_nc[i_cell]
        top_k_output = get_most_granular_top_k_calls(
            aggregated_scores, cl, min_acceptable_score, top_k, root_note, use_shortest_path
        )
        for k in range(top_k):
            top_k_calls_dict[f"{obs_prefix}_score_{k + 1}"].append(top_k_output[k][0])
            top_k_calls_dict[f"{obs_prefix}_name_{k + 1}"].append(top_k_output[k][2])
            top_k_calls_dict[f"{obs_prefix}_label_{k + 1}"].append(cl.cl_names_to_labels_map[top_k_output[k][2]])

    for k, v in top_k_calls_dict.items():
        adata.obs[k] = v


def get_knn_neighbor_indices(
    adata: AnnData,
    k_neighbors: int,
    representation_obsm_key: t.Optional[str] = None,
    compute_neighbors_if_missing: bool = True,
) -> t.List[np.ndarray]:
    """
    Return nearest query-cell neighbor indices for every cell.

    Uses the precomputed scanpy-style kNN graph stored in ``adata.obsp['distances']`` when present.
    If the graph is absent and ``compute_neighbors_if_missing`` is True, compute it with
    ``scanpy.pp.neighbors`` on ``adata.X`` (or the ``adata.obsm`` representation given by
    ``representation_obsm_key``). If the graph is absent and the flag is False, raise ``ValueError``.

    :param adata: AnnData object containing the query cells.
    :param k_neighbors: Number of nearest neighbors per cell to return.
    :param representation_obsm_key: Optional ``adata.obsm`` key holding the query-cell representation used to
        build the kNN graph when it is missing. ``None`` (default) uses ``adata.X``.
    :param compute_neighbors_if_missing: If True and no graph is stored in ``adata.obsp['distances']``,
        compute one with ``scanpy.pp.neighbors``. Defaults to True.

    :raises ValueError: If no graph exists and computation is disabled, or if the graph has fewer than
        ``k_neighbors`` non-self neighbors for any cell.

    :returns: A list with one array per cell, holding the nearest neighbor indices sorted by increasing distance.
    """
    if "distances" not in adata.obsp:
        if not compute_neighbors_if_missing:
            raise ValueError(
                "No precomputed kNN graph in adata.obsp['distances']. "
                "Set compute_neighbors_if_missing=True to compute one with scanpy.pp.neighbors."
            )
        logger.info("No precomputed kNN graph found; computing neighbors with scanpy.pp.neighbors.")
        import scanpy as sc

        # scanpy stores n_neighbors - 1 edges per cell (self is excluded), so request k_neighbors + 1
        # to end up with exactly k_neighbors neighbors per cell.
        sc.pp.neighbors(
            adata,
            n_neighbors=k_neighbors + 1,
            use_rep="X" if representation_obsm_key is None else representation_obsm_key,
        )
    else:
        logger.info("Found precomputed kNN graph in adata.obsp['distances']; using it.")

    distances_csr = adata.obsp["distances"].tocsr()
    neighbor_indices = []
    for i_cell in range(adata.n_obs):
        row = distances_csr[i_cell]
        order = np.argsort(row.data)
        # exclude the cell itself; scanpy graphs omit self, hand-built graphs may include it
        neighbor_idx = row.indices[order]
        neighbor_idx = neighbor_idx[neighbor_idx != i_cell]
        neighbor_indices.append(neighbor_idx[:k_neighbors])

    neighbor_counts = np.array([len(indices) for indices in neighbor_indices])
    if neighbor_counts.size and neighbor_counts.min() < k_neighbors:
        raise ValueError(
            f"Requested k_neighbors={k_neighbors}, but the precomputed kNN graph has only "
            f"{neighbor_counts.min()} to {neighbor_counts.max()} non-self neighbors per cell. "
            "Recompute the graph with more neighbors."
        )
    return neighbor_indices


def compute_most_granular_top_k_calls_knn(
    adata: AnnData,
    cl: CellOntologyCache,
    min_acceptable_score: float,
    k_neighbors: int = 5,
    representation_obsm_key: t.Optional[str] = None,
    compute_neighbors_if_missing: bool = True,
    aggregation_op: CellOntologyScoresAggregationOp = CellOntologyScoresAggregationOp.MEAN,
    aggregation_domain: CellOntologyScoresAggregationDomain = CellOntologyScoresAggregationDomain.ALL_CELLS,
    aggregation_score_threshold: float = 1e-4,
    top_k: int = 3,
    obs_prefix: str = "cas_knn_cell_type",
    root_note: str = CL_CELL_ROOT_NODE,
    use_shortest_path: bool = True,
):
    """
    Assign the most granular top-k cell type calls per kNN neighborhood in ``adata.obs``.

    Keeps single-cell resolution: for each cell, its CAS scores are aggregated over the cell itself
    plus its ``k_neighbors`` nearest query-cell neighbors (instead of over a hard Leiden cluster).
    Neighborhoods overlap, so every cell keeps its own refined annotation.

    :param adata: AnnData object with ``cas_cl_scores`` already inserted via :meth:`insert_ontology_aware_response`.
    :param cl: The CellOntologyCache object containing cell ontology term names and labels.
    :param min_acceptable_score: Minimum evidence score for a cell type call to be considered.
    :param representation_obsm_key: Optional ``adata.obsm`` key holding the query-cell representation used to
        build the kNN graph when it is missing. ``None`` (default) uses ``adata.X``.
    :param compute_neighbors_if_missing: If True and no kNN graph is stored in ``adata.obsp['distances']``,
        compute one with ``scanpy.pp.neighbors``. Defaults to True.
    :raises ValueError: If no graph exists and computation is disabled, or if the graph has fewer than
        ``k_neighbors`` non-self neighbors for any cell.
    :param aggregation_op: The aggregation operation to apply to the CAS scores within each neighborhood.
    :param aggregation_domain: The domain over which to perform the aggregation.
    :param aggregation_score_threshold: The threshold value for considering a CAS score as non-zero.
    :param top_k: Number of top calls to make per neighborhood.
    :param obs_prefix: Prefix for the ``.obs`` columns to write results into.
    :param root_note: Root node of the cell ontology used for ranking.
    :param use_shortest_path: Whether to use shortest (True) or longest (False) path depth for ranking.
    """
    neighbor_indices = get_knn_neighbor_indices(
        adata=adata,
        k_neighbors=k_neighbors,
        representation_obsm_key=representation_obsm_key,
        compute_neighbors_if_missing=compute_neighbors_if_missing,
    )

    n_obs = adata.n_obs
    top_k_calls_dict = dict()
    for k in range(top_k):
        top_k_calls_dict[f"{obs_prefix}_score_{k + 1}"] = [None] * n_obs
        top_k_calls_dict[f"{obs_prefix}_name_{k + 1}"] = [None] * n_obs
        top_k_calls_dict[f"{obs_prefix}_label_{k + 1}"] = [None] * n_obs

    for i_cell in range(n_obs):
        # self + its k nearest query-cell neighbors
        obs_indices = np.concatenate(([i_cell], neighbor_indices[i_cell]))
        aggregated_scores = get_aggregated_cas_ontology_aware_scores(
            adata, obs_indices, aggregation_op, aggregation_domain, aggregation_score_threshold
        )
        top_k_output = get_most_granular_top_k_calls(
            aggregated_scores, cl, min_acceptable_score, top_k, root_note, use_shortest_path
        )
        for k in range(top_k):
            top_k_calls_dict[f"{obs_prefix}_score_{k + 1}"][i_cell] = top_k_output[k][0]
            top_k_calls_dict[f"{obs_prefix}_name_{k + 1}"][i_cell] = top_k_output[k][2]
            top_k_calls_dict[f"{obs_prefix}_label_{k + 1}"][i_cell] = cl.cl_names_to_labels_map[top_k_output[k][2]]

    for k, v in top_k_calls_dict.items():
        adata.obs[k] = v


def compute_most_granular_top_k_calls_cluster(
    adata: AnnData,
    cl: CellOntologyCache,
    min_acceptable_score: float,
    cluster_label_obs_column: str,
    aggregation_op: CellOntologyScoresAggregationOp = CellOntologyScoresAggregationOp.MEAN,
    aggregation_domain: CellOntologyScoresAggregationDomain = CellOntologyScoresAggregationDomain.ALL_CELLS,
    aggregation_score_threshod: float = 1e-4,
    top_k: int = 3,
    obs_prefix: str = "cas_cell_type",
    root_note: str = CL_CELL_ROOT_NODE,
    use_shortest_path: bool = True,
):
    top_k_calls_dict = dict()
    for k in range(top_k):
        top_k_calls_dict[f"{obs_prefix}_score_{k + 1}"] = [None] * adata.n_obs
        top_k_calls_dict[f"{obs_prefix}_name_{k + 1}"] = [None] * adata.n_obs
        top_k_calls_dict[f"{obs_prefix}_label_{k + 1}"] = [None] * adata.n_obs

    def _update_list(target_list, indices, value):
        for idx in indices:
            target_list[idx] = value

    for cluster_label in adata.obs[cluster_label_obs_column].cat.categories:
        obs_indices = get_obs_indices_for_cluster(adata, cluster_label_obs_column, cluster_label)
        aggregated_scores = get_aggregated_cas_ontology_aware_scores(
            adata, obs_indices, aggregation_op, aggregation_domain, aggregation_score_threshod
        )
        top_k_output = get_most_granular_top_k_calls(
            aggregated_scores, cl, min_acceptable_score, top_k, root_note, use_shortest_path
        )
        for k in range(top_k):
            _update_list(top_k_calls_dict[f"{obs_prefix}_score_{k + 1}"], obs_indices, top_k_output[k][0])
            _update_list(top_k_calls_dict[f"{obs_prefix}_name_{k + 1}"], obs_indices, top_k_output[k][2])
            _update_list(
                top_k_calls_dict[f"{obs_prefix}_label_{k + 1}"],
                obs_indices,
                cl.cl_names_to_labels_map[top_k_output[k][2]],
            )

    for k, v in top_k_calls_dict.items():
        adata.obs[k] = v
