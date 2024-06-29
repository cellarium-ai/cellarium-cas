import typing as t

from cellarium.cas import exceptions

if t.TYPE_CHECKING:
    import anndata


def validate(
    adata: "anndata.AnnData",
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

    if adata_feature_schema_list == cas_feature_schema_list:
        return

    missing_features = len(cas_feature_schema_set - adata_feature_schema_set)
    extra_features = len(adata_feature_schema_set - cas_feature_schema_set)

    raise exceptions.DataValidationError(missing_features=missing_features, extra_features=extra_features)
