import typing as t

import numpy as np

from cellarium.cas import constants, exceptions

if t.TYPE_CHECKING:
    import anndata


def validate(
    adata: "anndata.AnnData",
    cas_feature_schema_list: t.List,
    feature_ids_column_name: str,
    count_matrix_input: constants.CountMatrixInput,
) -> None:
    """
    Validate input `anndata.AnnData` instance in concordance with feature schema `cas_feature_schema_list`
    and data is the correct type for compatibility with the backend code.
    Raise `exceptions.DataValidationError` with number of missing features and number of extra features
    that were present in dataset, along with a bool indicating whether the data is an incompatible type.

    :param adata: Instance to validate
    :param cas_feature_schema_list: List of features to be validated with
    :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. Default `index`.
    :param count_matrix_input: Where to obtain a feature expression count matrix from. Choice of: 'X', 'raw.X'
    """
    assert feature_ids_column_name == "index" or feature_ids_column_name in adata.var.columns.values, (
        "`feature_ids_column_name` should have a value of either 'index' "
        "or be present as a column in the `adata.var` object."
    )

    if feature_ids_column_name == "index":
        adata_feature_schema_list = adata.var.index.tolist()
    else:
        adata_feature_schema_list = adata.var[feature_ids_column_name].values.tolist()

    cas_feature_schema_set = set(cas_feature_schema_list)
    adata_feature_schema_set = set(adata_feature_schema_list)

    dtype_to_check = adata.X.dtype if count_matrix_input == constants.CountMatrixInput.X else adata.raw.X.dtype
    incompatible_x_type: t.Optional[str] = None if dtype_to_check == np.float32 else dtype_to_check

    incompatible_total_mrna_umis_type: t.Optional[str] = None
    if "total_mrna_umis" in adata.obs and adata.obs["total_mrna_umis"].dtype != np.float32:
        incompatible_total_mrna_umis_type = adata.obs["total_mrna_umis"].dtype

    if (
        adata_feature_schema_list == cas_feature_schema_list
        and incompatible_x_type is None
        and incompatible_total_mrna_umis_type is None
        # We want to be sure to sanitize if the count matrix is raw.X
        and count_matrix_input == constants.CountMatrixInput.X
    ):
        return

    missing_features = len(cas_feature_schema_set - adata_feature_schema_set)
    extra_features = len(adata_feature_schema_set - cas_feature_schema_set)

    raise exceptions.DataValidationError(
        missing_features=missing_features,
        extra_features=extra_features,
        incompatible_x_type=incompatible_x_type,
        incompatible_total_mrna_umis_type=incompatible_total_mrna_umis_type,
    )
