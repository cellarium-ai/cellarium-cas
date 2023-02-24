import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp

import warnings
from collections import defaultdict
from operator import itemgetter
from typing import Optional


def validate_adata_for_cas(
        adata: sc.AnnData,
        casp_feature_list_csv_path: str,
        int_count_matrix: str = 'X',
        gene_symbols_column_name: str = '__index__',
        gene_ids_column_name: str = 'gene_ids',
        missing_features_policy: str = 'replace_with_zero',
        extra_features_policy: str = 'ignore') -> sc.AnnData:
    
    # we only have the following policies implemented
    assert missing_features_policy == 'replace_with_zero'
    assert extra_features_policy == 'ignore'
    
    # where to find the integer count matrix?
    assert int_count_matrix in {'raw', 'X'}
    
    # fetch gene_symbols
    if gene_symbols_column_name == '__index__':
        adata_gene_symbols = adata.var.index.tolist()
    else:
        adata_gene_symbols = adata.var[gene_symbols_column_name].values.tolist()
    
    # fetch gene_ids
    if gene_ids_column_name == '__index__':
        adata_gene_ids = adata.var.index.tolist()
    else:
        adata_gene_ids = adata.var[gene_ids_column_name].values.tolist()
    
    # reformat anndata object for CAS
    casp_gene_id_list = pd.read_csv(casp_feature_list_csv_path, index_col=0)['gene_id'].values.tolist()
    adata_gene_id_list = adata_gene_ids
    gene_id_to_gene_symbol_map = {
        gene_id: gene_symbol for gene_symbol, gene_id in zip(adata_gene_symbols, adata_gene_ids)}
    casp_gene_symbol_list = [
        gene_id_to_gene_symbol_map[gene_id] if gene_id in gene_id_to_gene_symbol_map else 'N/A'
        for gene_id in casp_gene_id_list]
    
    gene_id_intersection = list(set(adata_gene_id_list).intersection(casp_gene_id_list))
    casp_gene_id_map = {gene_id: index for index, gene_id in enumerate(casp_gene_id_list)}
    adata_gene_id_map = {gene_id: index for index, gene_id in enumerate(adata_gene_id_list)}
    gene_id_intersection_casp_indices = list(map(casp_gene_id_map.get, gene_id_intersection))
    gene_id_intersection_adata_indices = list(map(adata_gene_id_map.get, gene_id_intersection))

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', sp.SparseEfficiencyWarning)
        n_cells = adata.shape[0]
        n_genes = len(casp_gene_id_list)
        X = sp.csc_matrix((n_cells, n_genes), dtype=np.float32)
        adata_X = adata.X if int_count_matrix == 'X' else adata.raw.X
        X[:, gene_id_intersection_casp_indices] = adata_X.tocsc()[:, gene_id_intersection_adata_indices]
        
        # workaround
        obs = adata.obs.copy()
        obs.index = pd.Index([str(x) for x in np.arange(adata.shape[0])])
        
        adata = sc.AnnData(
            X.tocsr(),
            obs=obs,
            obsm=adata.obsm,
            var=pd.DataFrame(data={'gene_symbols': casp_gene_symbol_list}, index=casp_gene_id_list))
    
    return adata


def reduce_cas_query_result_by_majority_vote(
        adata: sc.AnnData,
        cas_query_res: dict,
        output_cell_type_key: str = 'cas_cell_type',
        output_cell_type_confidence_score_key: str = 'cas_cell_type_confidence_score'):

    majority_vote_cell_type_list = []
    majority_vote_confidence_score_list = []
    for single_cell_query in cas_query_res:
        total_cell_count = 0
        best_cell_count = -1
        best_cell_type = None
        for match in single_cell_query['matches']:
            total_cell_count += match['cell_count']
            if match['cell_count'] > best_cell_count:
                best_cell_count = match['cell_count']
                best_cell_type = match['cell_type']
        majority_vote_cell_type_list.append(best_cell_type)
        majority_vote_confidence_score_list.append(best_cell_count / total_cell_count)

    adata.obs[output_cell_type_key] = majority_vote_cell_type_list
    adata.obs[output_cell_type_confidence_score_key] = majority_vote_confidence_score_list


def reduce_cas_query_result_by_min_distance(
        adata: sc.AnnData,
        cas_query_res: dict,
        output_cell_type_key: str = 'cas_cell_type',
        output_cell_type_min_distance_key: str = 'cas_cell_type_min_distance'):
    
    min_distance_cell_type_list = []
    min_distance_list = []
    for single_cell_query in cas_query_res:
        min_distance = np.inf
        best_cell_type = None
        for match in single_cell_query['matches']:
            if match['min_distance'] < min_distance:
                min_distance = match['min_distance']
                best_cell_type = match['cell_type']
        min_distance_cell_type_list.append(best_cell_type)
        min_distance_list.append(min_distance)
        
    adata.obs[output_cell_type_key] = min_distance_cell_type_list
    adata.obs[output_cell_type_min_distance_key] = min_distance_list

    
def reduce_cas_query_result_by_majority_vote_per_cluster(
        adata: sc.AnnData,
        cas_query_res: dict,
        cluster_key: str = 'leiden',
        output_cell_type_key: str = 'cas_per_cluster_cell_type',
        output_cell_type_confidence_score_key: str = 'cas_per_cluster_cell_type_confidence_score',
        ignore_set: set = set()) -> dict:

    assert adata.shape[0] == len(cas_query_res)
    assert cluster_key in adata.obs

    output_cell_type_list = [None] * adata.shape[0]
    output_cell_type_confidence_score_list = [None] * adata.shape[0]
    cluster_detailed_info_dict = dict()
    
    for cluster_id in adata.obs[cluster_key].values.categories:

        # obtain cell indices belonging to cluster_id
        cluster_cell_indices = np.where(adata.obs[cluster_key] == cluster_id)[0]
        assert len(cluster_cell_indices) > 0

        # summarize the hits of all cells
        cluster_query = defaultdict(int)
        for cell_index in cluster_cell_indices:
            single_cell_query = cas_query_res[cell_index]
            for match in single_cell_query['matches']:
                if match['cell_type'] in ignore_set:
                    continue
                cluster_query[match['cell_type']] += match['cell_count']

        # identify best cell type for the cluster
        total_cell_count = sum(cluster_query.values())
        sorted_hits_and_freqs = sorted(
            [(cell_type, count / total_cell_count) for cell_type, count in cluster_query.items()],
            key=itemgetter(1),
            reverse=True)
        best_cell_type = sorted_hits_and_freqs[0][0]
        best_cell_type_confidence_score = sorted_hits_and_freqs[0][1]

        # bookkeeping
        cluster_detailed_info_dict[cluster_id] = sorted_hits_and_freqs        
        for cell_index in cluster_cell_indices:
            output_cell_type_list[cell_index] = best_cell_type
            output_cell_type_confidence_score_list[cell_index] = best_cell_type_confidence_score
            
    adata.obs[output_cell_type_key] = output_cell_type_list
    adata.obs[output_cell_type_confidence_score_key] = output_cell_type_confidence_score_list
    
    return cluster_detailed_info_dict
    
    