import anndata
import numpy as np
import pandas as pd

from cellarium.cas.data_preparation import callbacks


def test_calculate_total_mrna_umis():
    """
    Test :func:`data_preparation.callbacks.calculate_total_mrna_umis` function.

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

    result_adata = callbacks.calculate_total_mrna_umis(adata)

    assert (
        callbacks.TOTAL_MRNA_UMIS_COLUMN_NAME in result_adata.obs.columns
    ), f"'{callbacks.TOTAL_MRNA_UMIS_COLUMN_NAME}' column not found in .obs"

    np.testing.assert_array_equal(
        result_adata.obs[callbacks.TOTAL_MRNA_UMIS_COLUMN_NAME],
        [3, 1, 3],
        err_msg="Total mRNA UMI calculations are incorrect",
    )

    # Check if the original AnnData object is unmodified (immutability test)
    assert adata.shape == result_adata.shape, "Original AnnData object was modified"
    assert "total_mrna_umis" not in adata.obs.columns, "Original AnnData object was modified"
