import typing as t

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

    if (
        adata_feature_schema_list == cas_feature_schema_list
        # We want to be sure to sanitize if the count matrix is raw.X
        and count_matrix_input == constants.CountMatrixInput.X
    ):
        return

    missing_features = len(cas_feature_schema_set - adata_feature_schema_set)
    extra_features = len(adata_feature_schema_set - cas_feature_schema_set)

    raise exceptions.DataValidationError(
        missing_features=missing_features,
        extra_features=extra_features,
    )
