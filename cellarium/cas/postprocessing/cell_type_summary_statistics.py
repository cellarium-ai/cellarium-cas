from __future__ import annotations

import typing as t
from collections import defaultdict
from operator import itemgetter

import numpy as np
from anndata import AnnData
from tqdm.notebook import tqdm

ALLOWED_WNN_STRATEGIES = {"connectivities"}


def reduce_cas_query_result_by_majority_vote(
    adata: AnnData,
    cas_cell_type_summary_statistics_response: list,
    output_cell_type_obs_column: str = "cas_cell_type",
    output_cell_type_confidence_score_obs_column: str = "cas_cell_type_confidence_score",
):
    """Reduce the CAS cell type summary statistics response by majority vote.

    This function takes the result of a CAS cell type summary statistics response (the output of
    :class:`cellarium.cas.client.annotate_matrix_cell_type_summary_statistics_strategy`) and reduces
    it by performing a majority vote to determine the most likely cell type for each query. It also
    calculates a confidence score based on the ratio of the best cell count to the total cell count.

    :param adata: An Annotated Data object containing the query result.
    :type adata: anndata.AnnData
    :param cas_cell_type_summary_statistics_response: The summary statistics response from the CAS query.
    :type cas_cell_type_summary_statistics_response: list
    :param output_cell_type_obs_column: The name of the observation column to store the predicted cell types.
    :type output_cell_type_obs_column: str, optional
    :param output_cell_type_confidence_score_obs_column: The name of the observation column to store the
      confidence scores.
    :type output_cell_type_confidence_score_obs_column: str, optional
    """

    majority_vote_cell_type_list = []
    majority_vote_confidence_score_list = []
    for single_cell_query in cas_cell_type_summary_statistics_response:
        total_cell_count = 0
        best_cell_count = -1
        best_cell_type = None
        for match in single_cell_query["matches"]:
            total_cell_count += match["cell_count"]
            if match["cell_count"] > best_cell_count:
                best_cell_count = match["cell_count"]
                best_cell_type = match["cell_type"]
        majority_vote_cell_type_list.append(best_cell_type)
        majority_vote_confidence_score_list.append(best_cell_count / total_cell_count)

    adata.obs[output_cell_type_obs_column] = majority_vote_cell_type_list
    adata.obs[output_cell_type_confidence_score_obs_column] = majority_vote_confidence_score_list


def reduce_cas_cell_type_summary_statistics_response_by_min_distance(
    adata: AnnData,
    cas_cell_type_summary_statistics_response: list,
    output_cell_type_obs_column: str = "cas_cell_type",
    output_cell_type_min_distance_obs_column: str = "cas_cell_type_min_distance",
):
    """
    Reduce the CAS query result by selecting the cell type with the minimum distance for each query.

    :param adata: Annotated data object.
    :type adata: AnnData
    :param cas_cell_type_summary_statistics_response: CAS response.
    :type cas_cell_type_summary_statistics_response: list
    :param output_cell_type_obs_column: Name of the output observation column for cell types.
      (default: "cas_cell_type")
    :type output_cell_type_obs_column: str
    :param output_cell_type_min_distance_obs_column: Name of the output observation column for minimum distances.
      (default: "cas_cell_type_min_distance")
    :type output_cell_type_min_distance_obs_column: str
    :return: None
    """

    min_distance_cell_type_list = []
    min_distance_list = []
    for single_cell_query in cas_cell_type_summary_statistics_response:
        min_distance = np.inf
        best_cell_type = None
        for match in single_cell_query["matches"]:
            if match["min_distance"] < min_distance:
                min_distance = match["min_distance"]
                best_cell_type = match["cell_type"]
        min_distance_cell_type_list.append(best_cell_type)
        min_distance_list.append(min_distance)

    adata.obs[output_cell_type_obs_column] = min_distance_cell_type_list
    adata.obs[output_cell_type_min_distance_obs_column] = min_distance_list


def reduce_cas_cell_type_summary_statistics_response_by_majority_vote_per_cluster(
    adata: AnnData,
    cas_cell_type_summary_statistics_response: dict,
    cluster_key: str = "leiden",
    output_cell_type_obs_column: str = "cas_per_cluster_cell_type",
    output_cell_type_confidence_score_obs_column: str = "cas_per_cluster_cell_type_confidence_score",
    ignore_set: t.Optional[set] = None,
) -> dict:
    """
    Reduce the CAS query result by performing a majority vote per cluster.

    This function takes the result of a CAS query and reduces it by performing a majority vote per cluster.
    It calculates the most likely cell type for each cluster based on the majority vote and assigns it to the
    corresponding cells in the query result. It also calculates a confidence score based on the ratio of the
    best cell count to the total cell count.

    :param adata: An Annotated Data object containing the query result.
    :type adata: anndata.AnnData
    :param cas_query_res: The CAS query result dictionary.
    :type cas_query_res: dict
    :param cluster_key: The key for the cluster column in `adata.obs`.
    :type cluster_key: str, optional
    :param output_cell_type_obs_column: The name of the observation column to store the predicted cell types.
    :type output_cell_type_obs_column: str, optional
    :param output_cell_type_confidence_score_obs_column: The name of the observation column to store the
      confidence scores.
    :type output_cell_type_confidence_score_obs_column: str, optional
    :param ignore_set: A set of cell types to ignore during the majority vote.
    :type ignore_set: set, optional
    :return: A dictionary containing detailed information about the majority vote for each cluster.
    :rtype: dict
    """
    if len(adata) != len(cas_cell_type_summary_statistics_response):
        raise ValueError("Lengths of `cas_query_res` and `adata` should be the same")
    if cluster_key not in adata.obs:
        raise ValueError("`cluster_key` does not correspond to a column in `adata.obs`")
    if ignore_set is None:
        ignore_set = set()
    output_cell_type_list = [None] * adata.shape[0]
    output_cell_type_confidence_score_list = [None] * adata.shape[0]
    cluster_detailed_info_dict = dict()
    for cluster_id in adata.obs[cluster_key].values.categories:
        cluster_cell_indices = np.where(adata.obs[cluster_key] == cluster_id)[0]
        if len(cluster_cell_indices) == 0:
            raise ValueError(f"No cells found belonging to cluster_id {cluster_id}")
        cluster_query = defaultdict(int)
        for cell_index in cluster_cell_indices:
            single_cell_query = cas_cell_type_summary_statistics_response[cell_index]
            for match in single_cell_query["matches"]:
                if match["cell_type"] in ignore_set:
                    continue
                cluster_query[match["cell_type"]] += match["cell_count"]
        total_cell_count = sum(cluster_query.values())
        sorted_hits_and_freqs = sorted(
            [(cell_type, count / total_cell_count) for cell_type, count in cluster_query.items()],
            key=itemgetter(1),
            reverse=True,
        )
        best_cell_type = sorted_hits_and_freqs[0][0]
        best_cell_type_confidence_score = sorted_hits_and_freqs[0][1]
        cluster_detailed_info_dict[cluster_id] = sorted_hits_and_freqs
        for cell_index in cluster_cell_indices:
            output_cell_type_list[cell_index] = best_cell_type
            output_cell_type_confidence_score_list[cell_index] = best_cell_type_confidence_score
    adata.obs[output_cell_type_obs_column] = output_cell_type_list
    adata.obs[output_cell_type_confidence_score_obs_column] = output_cell_type_confidence_score_list

    return cluster_detailed_info_dict


def _get_weights_via_fuzzy_simplicial_sets(
    i_obs: int,
    adata: AnnData,
    n_neighbors: int,
    self_connectivity: float,
    connectivities_key: str,
) -> tuple:
    assert connectivities_key in adata.obsp, "Please run `scanpy.pp.neighbors` first"
    connectivity_values = adata.obsp[connectivities_key][i_obs].data
    connectivity_indices = adata.obsp[connectivities_key][i_obs].indices

    # append self
    connectivity_values = np.append(connectivity_values, self_connectivity)
    connectivity_indices = np.append(connectivity_indices, i_obs)

    # convert to weights
    _order = np.argsort(connectivity_values)[::-1][:n_neighbors]
    sorted_connectivity_values = connectivity_values[_order]
    sorted_connectivity_indices = connectivity_indices[_order]
    weights = sorted_connectivity_values / np.sum(sorted_connectivity_values)
    return sorted_connectivity_indices, weights


def _get_cell_type_probs(single_cell_query: t.Dict, all_cell_types: t.List, cell_type_to_idx_map: t.Dict) -> np.ndarray:
    probs = np.zeros((len(all_cell_types),))
    for match in single_cell_query["matches"]:
        probs[cell_type_to_idx_map[match["cell_type"]]] += match["cell_count"]
    return probs / np.sum(probs)


def reduce_cas_cell_type_summary_statistics_response_by_wnn(
    adata: AnnData,
    cas_cell_type_summary_statistics_response: dict,
    n_neighbors: int = 10,
    wnn_strategy: str = "connectivities",
    connectivities_key: str = "connectivities",
    self_connectivity: float = 1.0,
    min_n_cells_per_type: int = 10,
    output_unreliable_type: str = "Unknown or Unconfident",
    output_cell_type_obs_column: str = "cas_cell_type",
    output_cell_type_confidence_score_obs_column: str = "cas_cell_type_confidence_score",
):
    """
    Reduce cell type summary statistics response by Weighted Nearest Neighbors (WNN).

    Args:
        adata (AnnData): Annotated data object.
        cas_cell_type_summary_statistics_response (dict): Cell type summary statistics response from CAS.
        n_neighbors (int, optional): Number of neighbors to consider. Defaults to 10.
        wnn_strategy (str, optional): WNN strategy. Defaults to "connectivities".
        connectivities_key (str, optional): Key for connectivities in `adata.uns['neighbors']`.
          Defaults to "connectivities".
        self_connectivity (float, optional): Self-connectivity value. Defaults to 1.0.
        min_n_cells_per_type (int, optional): Minimum number of cells per cell type.
          Defaults to 10.
        output_unreliable_type (str, optional): Output value for unreliable cell types.
          Defaults to "Unknown or Unconfident".
        output_cell_type_obs_column (str, optional): Output column name for cell type. Defaults to "cas_cell_type".
        output_cell_type_confidence_score_obs_column (str, optional): Output column name for cell type confidence score. Defaults to "cas_cell_type_confidence_score".

    Returns:
        tuple: A tuple containing the cell type to index mapping and a list of cell type probability arrays.
    """

    if wnn_strategy not in ALLOWED_WNN_STRATEGIES:
        raise ValueError(
            f"`wnn_strategy` should be one of {', '.join(ALLOWED_WNN_STRATEGIES)}, got {wnn_strategy} instead"
        )
    if n_neighbors > adata.uns["neighbors"]["params"]["n_neighbors"]:
        raise ValueError(
            "`n_neighbors` should be less than or equal to `adata.uns['neighbors']['params']['n_neighbors']`"
        )
    if len(adata) != len(cas_cell_type_summary_statistics_response):
        raise ValueError("Lengths of `cas_query_res` and `adata` should be the same")

    all_cell_types = []
    for single_cell_query in cas_cell_type_summary_statistics_response:
        for match in single_cell_query["matches"]:
            all_cell_types.append(match["cell_type"])
    all_cell_types = list(set(all_cell_types))
    cell_type_to_idx_map = {cell_type: idx for idx, cell_type in enumerate(all_cell_types)}

    cell_type_probs_list = []
    majority_vote_cell_type_list = []
    majority_vote_confidence_score_list = []

    n_cells = len(cas_cell_type_summary_statistics_response)
    for i in tqdm(range(n_cells)):
        neighbor_indices, neighbor_weights = _get_weights_via_fuzzy_simplicial_sets(
            i,
            adata,
            n_neighbors,
            self_connectivity,
            connectivities_key=connectivities_key,
        )
        if not np.isclose(np.sum(neighbor_weights), 1.0):
            raise ValueError("Sum of Neighbor weights is not close to 1")

        cell_type_probs = np.zeros((len(all_cell_types),))
        for j, weight in zip(neighbor_indices, neighbor_weights):
            cell_type_probs += weight * _get_cell_type_probs(
                cas_cell_type_summary_statistics_response[j], all_cell_types, cell_type_to_idx_map
            )

        best_cell_type_idx = np.argmax(cell_type_probs)
        best_cell_type = all_cell_types[best_cell_type_idx]
        best_cell_type_prob = cell_type_probs[best_cell_type_idx]

        cell_type_probs_list.append(cell_type_probs)
        majority_vote_cell_type_list.append(best_cell_type)
        majority_vote_confidence_score_list.append(best_cell_type_prob)

    all_cell_types = set(majority_vote_cell_type_list)
    update_map = dict()
    for cell_type in all_cell_types:
        n_counts = sum([_cell_type == cell_type for _cell_type in majority_vote_cell_type_list])
        if n_counts < min_n_cells_per_type:
            update_map[cell_type] = output_unreliable_type
        else:
            update_map[cell_type] = cell_type
    majority_vote_cell_type_list = list(map(update_map.get, majority_vote_cell_type_list))

    adata.obs[output_cell_type_obs_column] = majority_vote_cell_type_list
    adata.obs[output_cell_type_confidence_score_obs_column] = majority_vote_confidence_score_list

    return cell_type_to_idx_map, cell_type_probs_list


def hex_to_rgb(value: str) -> np.ndarray:
    """
    Convert a hexadecimal color value to an RGB array.

    :param value: The hexadecimal color value to convert.
    :type value: str
    :return: The RGB array representing the color.
    :rtype: np.ndarray
    """
    value = value.lstrip("#")
    lv = len(value)
    return np.asarray(tuple(float(int(value[i : i + lv // 3], 16)) for i in range(0, lv, lv // 3)))


def rgb_to_tuple(value: np.ndarray) -> tuple:
    """
    Convert an RGB value represented as a NumPy array to a tuple of integers.

    :param value: The RGB value as a NumPy array.
    :type value: np.ndarray
    :return: The RGB value as a tuple of integers.
    :rtype: tuple
    """
    return int(value[0]), int(value[1]), int(value[2])


def rgb_to_hex(r: int, g: int, b: int) -> str:
    """
    Convert RGB values to a hexadecimal color code.

    :param r: The red component of the RGB color (0-255).
    :param g: The green component of the RGB color (0-255).
    :param b: The blue component of the RGB color (0-255).
    :return: The hexadecimal color code representing the RGB color.
    """
    return "#%02x%02x%02x" % (r, g, b)


def get_interpolated_cell_type_colors(
    adata: AnnData,
    cell_type_to_idx_map: t.Dict[str, int],
    cell_type_probs_list: t.List[np.ndarray],
    na_cell_type_key: str = "Unknown or Unconfident",
) -> np.ndarray:
    """Get interpolated cell type colors based on the cell type probabilities.

    This function takes the cell type probabilities and maps them to the cell type colors. The cell type colors are
    interpolated based on the probabilities.

    :param adata: An Annotated Data object containing the query result.
    :type adata: anndata.AnnData
    :param cell_type_to_idx_map: A dictionary mapping cell types to indices.
    :type cell_type_to_idx_map: dict
    :param cell_type_probs_list: A list of cell type probabilities.
    :type cell_type_probs_list: list
    :param na_cell_type_key: The key for the unknown cell type.
    :type na_cell_type_key: str, optional
    :return: The interpolated cell type colors.
    :rtype: np.ndarray
    """
    idx_to_cell_type_map: t.Dict[int, str] = {idx: cell_type for cell_type, idx in cell_type_to_idx_map.items()}
    all_cell_types: t.List[str] = list(map(idx_to_cell_type_map.get, range(len(cell_type_to_idx_map))))
    plot_cell_type_colors: np.ndarray = adata.uns["cas_cell_type_colors"]
    plot_cell_types: t.List[str] = list(adata.obs["cas_cell_type"].values.categories)
    plot_cell_types_to_idx_map: t.Dict[str, int] = {cell_type: idx for idx, cell_type in enumerate(plot_cell_types)}

    cell_type_probs_nk: np.ndarray = np.asarray(cell_type_probs_list)
    cell_type_map_kq: np.ndarray = np.zeros((len(all_cell_types), len(plot_cell_types)))
    na_index: int = plot_cell_types.index(na_cell_type_key)
    plot_cell_type_colors[na_index] = rgb_to_hex(r=255, g=255, b=255)
    for k, input_cell_type in enumerate(all_cell_types):
        if input_cell_type in plot_cell_types_to_idx_map:
            q: int = plot_cell_types_to_idx_map[input_cell_type]
        else:
            q: int = na_index
        cell_type_map_kq[k, q] += 1

    plot_colors_q3: np.ndarray = np.asarray(list(hex_to_rgb(hex_color) for hex_color in plot_cell_type_colors))
    cell_colors_n3: np.ndarray = (cell_type_probs_nk @ cell_type_map_kq @ plot_colors_q3) / 255

    return np.clip(cell_colors_n3, 0.0, 1.0)
