import anndata
import numpy as np
import pandas as pd

from cellarium.cas import constants
from cellarium.cas.preprocessing import callbacks

np_random_state = np.random.RandomState(0)


def test_calculate_total_mrna_umis_X():
    """
    Test :func:`preprocessing.callbacks.calculate_total_mrna_umis` function.

    Assert the function indeed creates a ``"total_mrna_umis"`` obs column, check the values correspond to the sum
    over the axis, and insure the input :class:`anndata.AnnData` instance remains unchanged.

    Raises:
    - AssertionError: If either the output :class:`anndata.AnnData` instance doesn't have correct numbers in
        ``.obs["total_mrna_umis"]`` or if the input :class:`anndata.AnnData` instance has changed.
    """
    X = np.array(
        [
            [1, 0, 2],  # Cell 1: total mRNA UMIs = 3
            [0, 0, 1],  # Cell 2: total mRNA UMIs = 1
            [2, 1, 0],  # Cell 3: total mRNA UMIs = 3
        ]
    )

    obs = pd.DataFrame(index=["cell1", "cell2", "cell3"])
    adata = anndata.AnnData(X=X, obs=obs)

    callbacks.calculate_total_mrna_umis(adata, count_matrix_input=constants.CountMatrixInput.X)

    assert (
        callbacks.TOTAL_MRNA_UMIS_COLUMN_NAME in adata.obs.columns
    ), f"'{callbacks.TOTAL_MRNA_UMIS_COLUMN_NAME}' column not found in .obs"

    np.testing.assert_array_equal(
        adata.obs[callbacks.TOTAL_MRNA_UMIS_COLUMN_NAME],
        [3, 1, 3],
        err_msg="Total mRNA UMI calculations are incorrect",
    )


def test_calculate_total_mrna_umis_raw_X():
    """
    Test :func:`preprocessing.callbacks.calculate_total_mrna_umis` function when X has normalized data
    and ``count_matrix_name="raw.X"``.

    Assert the function indeed creates a ``"total_mrna_umis"`` obs column, check the values correspond to the sum
    over the axis in raw count matrix, and insure the input :class:`anndata.AnnData` instance remains unchanged.

    Raises:
    - AssertionError: If either the output :class:`anndata.AnnData` instance doesn't have correct numbers in
        ``.obs["total_mrna_umis"]`` or if the input :class:`anndata.AnnData` instance has changed.
    """
    X = np.array(
        [
            [1, 0, 2],  # Cell 1: total mRNA UMIs = 3
            [0, 0, 1],  # Cell 2: total mRNA UMIs = 1
            [2, 1, 0],  # Cell 3: total mRNA UMIs = 3
        ]
    )

    obs = pd.DataFrame(index=["cell1", "cell2", "cell3"])
    adata = anndata.AnnData(X=X, obs=obs)
    # We need to have an anndata file that has normalized values in `X` and raw values in `raw.X`
    adata.raw = adata
    adata.X = 10_000 * adata.X / (adata.X.sum(axis=1) + 0.000001)  # Using Normalize Total for matrix X

    callbacks.calculate_total_mrna_umis(adata, count_matrix_input=constants.CountMatrixInput.RAW_X)

    assert (
        callbacks.TOTAL_MRNA_UMIS_COLUMN_NAME in adata.obs.columns
    ), f"'{callbacks.TOTAL_MRNA_UMIS_COLUMN_NAME}' column not found in .obs"

    np.testing.assert_array_equal(
        adata.obs[callbacks.TOTAL_MRNA_UMIS_COLUMN_NAME],
        [3, 1, 3],
        err_msg="Total mRNA UMI calculations are incorrect",
    )


def test_ensure_matrix_is_float32_X():
    """
    Test :func:`preprocessing.callbacks.ensure_matrix_is_float32` function when X is not of type float32.
    """
    X = np_random_state.randn(5, 5).astype(np.float64)
    obs = pd.DataFrame(index=["cell1", "cell2", "cell3", "cell4", "cell5"])
    adata = anndata.AnnData(X=X, obs=obs)

    callbacks.ensure_matrix_is_float32(adata, count_matrix_input=constants.CountMatrixInput.X)

    assert adata.X.dtype == np.float32, "X matrix is not of type float32"


def test_ensure_matrix_is_float32_raw_X():
    """
    Test :func:`preprocessing.callbacks.ensure_matrix_is_float32` function when X is not of type float32.
    """
    raw_X = np_random_state.randn(5, 5).astype(np.int32)
    obs = pd.DataFrame(index=["cell1", "cell2", "cell3", "cell4", "cell5"])
    adata = anndata.AnnData(X=raw_X, obs=obs)
    adata.raw = adata

    callbacks.ensure_matrix_is_float32(adata, count_matrix_input=constants.CountMatrixInput.RAW_X)

    assert adata.raw.X.dtype == np.float32, "raw X matrix is not of type float32"
