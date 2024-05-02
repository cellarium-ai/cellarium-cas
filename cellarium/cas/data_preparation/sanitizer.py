import typing as t

import anndata
import numpy as np
import pandas as pd
import scipy.sparse as sp

from cellarium.cas import constants, exceptions
from cellarium.cas.data_preparation import callbacks


def _get_adata_var_index_or_by_column(adata: anndata.AnnData, var_column_name: str) -> t.List[str]:
    if var_column_name == "index":
        return adata.var.index.tolist()

    return adata.var[var_column_name].values.tolist()


def validate(
    adata: anndata.AnnData,
    cas_feature_schema_list: t.List,
    feature_ids_column_name: str,
) -> None:
    """
    Validate input `anndata.AnnData` instance in concordance with feature schema `cas_feature_schema_list`.
    Raise `exceptions.DataValidationError` with number of missing features and number of extra features
    that were present in dataset.

    :param adata: Instance to validate
    :param cas_feature_schema_list: List of features to be validated with
    :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. Default `index`.
    """
    if feature_ids_column_name != "index" and feature_ids_column_name not in adata.var.columns.values:
        raise ValueError(
            "`feature_ids_column_name` should have a value of either 'index' "
            "or be present as a column in the `adata.var` object."
        )

    adata_feature_schema_list = _get_adata_var_index_or_by_column(adata=adata, var_column_name=feature_ids_column_name)

    cas_feature_schema_set = set(cas_feature_schema_list)
    adata_feature_schema_set = set(adata_feature_schema_list)

    if adata_feature_schema_list == cas_feature_schema_list:
        return

    missing_features = len(cas_feature_schema_set - adata_feature_schema_set)
    extra_features = len(adata_feature_schema_set - cas_feature_schema_set)

    raise exceptions.DataValidationError(missing_features=missing_features, extra_features=extra_features)


def sanitize(
    adata: anndata.AnnData,
    cas_feature_schema_list: t.List[str],
    count_matrix_input: constants.CountMatrixInput,
    feature_ids_column_name: str,
    feature_names_column_name: t.Optional[str] = None,
) -> anndata.AnnData:
    """
    Cellarium CAS sanitizing script. Returns a new `anndata.AnnData` instance, based on the feature expression
    matrix of the input instance. Extra features get omitted. Missing features get filled with zeros.

    :param adata: Instance to sanitize
    :param cas_feature_schema_list: List of Ensembl feature ids to rely on.
    :param count_matrix_input: Where to obtain a feature expression count matrix from. Choice of: 'X', 'raw.X'
    :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. Default `index`.
    :param feature_names_column_name: Column name where to obtain feature names. If not provided, no feature names
        should be mapped |br|
    `Default`: ``None``
    :return: `anndata.AnnData` instance that is ready to use with CAS
    """
    if feature_ids_column_name != "index" and feature_ids_column_name not in adata.var.columns.values:
        raise ValueError(
            "`feature_ids_column_name` should have a value of either 'index' "
            "or be present as a column in the `adata.var` object."
        )
    if feature_names_column_name not in {"index", None} and feature_names_column_name not in adata.var.columns.values:
        raise ValueError(
            "`feature_ids_name_column_name` should have a value of either 'index' "
            "or be present as a column in the `adata.var` object."
        )

    callbacks.pre_sanitize_callback(adata=adata, count_matrix_input=count_matrix_input)

    adata_feature_schema_list = _get_adata_var_index_or_by_column(adata=adata, var_column_name=feature_ids_column_name)
    original_obs_ids = adata.obs.index.values

    cas_feature_schema_set = set(cas_feature_schema_list)
    adata_feature_schema_set = set(adata_feature_schema_list)
    feature_id_intersection = adata_feature_schema_set.intersection(cas_feature_schema_set)

    cas_feature_id_map = {feature_id: index for index, feature_id in enumerate(cas_feature_schema_list)}
    adata_feature_id_map = {feature_id: index for index, feature_id in enumerate(adata_feature_schema_list)}
    feature_id_intersection_cas_indices = list(map(cas_feature_id_map.get, feature_id_intersection))
    feature_id_intersection_adata_indices = list(map(adata_feature_id_map.get, feature_id_intersection))

    n_cells = adata.shape[0]
    n_features = len(cas_feature_schema_list)
    input_matrix = adata.X if count_matrix_input == constants.CountMatrixInput.X else adata.raw.X

    # Translate the columns from one matrix to another, convert to COO format to make this efficient.
    col_trans = np.zeros(n_features, dtype=int)
    for i, k in enumerate(feature_id_intersection_cas_indices):
        col_trans[i] = k
    vals = input_matrix.tocsc()[:, feature_id_intersection_adata_indices]
    vals = vals.tocoo()
    new_col = col_trans[vals.col]
    result_matrix = sp.coo_matrix((vals.data, (vals.row, new_col)), shape=(n_cells, n_features))
    del col_trans, vals, new_col

    # Create `obs` index
    obs = adata.obs.copy()
    obs.index = original_obs_ids

    var_df_data = {}
    if feature_names_column_name is not None:
        adata_feature_name_list = _get_adata_var_index_or_by_column(
            adata=adata, var_column_name=feature_names_column_name
        )
        gene_id_to_gene_symbol_map = {
            gene_id: gene_symbol for gene_symbol, gene_id in zip(adata_feature_name_list, adata_feature_schema_list)
        }

        cas_feature_name_list = [
            gene_id_to_gene_symbol_map[gene_id] if gene_id in gene_id_to_gene_symbol_map else "N/A"
            for gene_id in cas_feature_schema_list
        ]

        var_df_data["feature_name"] = cas_feature_name_list

    var_df = pd.DataFrame(index=cas_feature_schema_list, data=var_df_data)
    return anndata.AnnData(
        result_matrix.tocsr(),
        obs=obs,
        obsm=adata.obsm,
        var=var_df,
    )
