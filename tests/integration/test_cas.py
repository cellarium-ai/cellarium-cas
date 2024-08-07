"""
Test Cell Annotation (Full Cycle)

The purpose of this case is to test whether the CAS client tool functions properly. It connects to the backend,
and annotates a test dataset using the default model.
"""

import anndata
import numpy as np

from cellarium.cas import CASClient
from tests.integration import constants

np_random_state = np.random.RandomState(0)


def test_cell_annotation(test_api_token: str, test_api_url: str):
    """
    Test the cell annotation functionality of the CASClient by reading a sample dataset and using the
    annotate_anndata method. It verifies that the number of annotations matches the original cell count.

    Parameters:
    - test_api_token (str): The API token used to authenticate with the CASClient.
    - test_api_url (str): The API url of the CAS backend that the CASClient connects to.

    Raises:
    - AssertionError: If the number of annotations doesn't match the original cell count in the dataset.

    Notes:
    - Ensure that TEST_ADATA_PATH is correctly defined and points to a valid dataset file.
    """
    cas = CASClient(api_token=test_api_token, api_url=test_api_url)
    adata = anndata.read_h5ad(constants.TEST_ADATA_PATH)
    result = cas.annotate_anndata(adata=adata, chunk_size=50)
    assert len(result) == len(adata), "Result length does not correspond to original number of cells"
    # Make sure that the results were all successful
    error_results = []
    for i, r in enumerate(result):
        if r["matches"] == []:
            error_results.append(i)
    assert len(error_results) == 0, f"Error in cell annotation at indices {error_results}"
