from typing import Any

import numpy as np
import pandas as pd
from anndata import AnnData


def get_obs_indices_for_cluster(adata: AnnData, obs_column: str, value: Any) -> np.ndarray:
    """
    Get the indices of observations in a specific cluster.

    :param adata: Annotated data object.
    :type adata: AnnData
    :param obs_column: Name of the observation column to filter on.
    :type obs_column: str
    :param value: Value to filter the observation column on.
    :type value: Any
    :return: Array of indices of observations in the specified cluster.
    :rtype: np.ndarray
    """
    assert obs_column in adata.obs.columns
    assert isinstance(adata.obs[obs_column].dtype, pd.CategoricalDtype)
    assert value in adata.obs[obs_column].cat.categories
    return (adata.obs[obs_column] == value).values.nonzero()[0]
