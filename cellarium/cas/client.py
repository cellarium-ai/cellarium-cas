import asyncio
import datetime
import functools
import math
import operator
import time
import typing as t

import anndata

from cellarium.cas import _io, data_preparation, exceptions, service

NUM_ATTEMPTS_PER_CHUNK_DEFAULT = 3


class CASClient:
    """
    Service that is designed to communicate with the Cellarium Cloud Backend.

    :param api_token: API token issued by the Cellarium team
    :param num_attempts_per_chunk: Number of attempts the client should make to annotate each chunk. |br|
        `Default:` ``3``
    """

    def __init__(self, api_token: str, num_attempts_per_chunk: int = NUM_ATTEMPTS_PER_CHUNK_DEFAULT) -> None:
        self._print("Connecting to the Cellarium Cloud backend...")
        self.cas_api_service = service.CASAPIService(api_token=api_token)
        self.feature_schemas = []
        self._feature_schemas_cache = {}
        self.num_attempts_per_chunk = num_attempts_per_chunk

    @staticmethod
    def _get_number_of_chunks(adata, chunk_size):
        return math.ceil(len(adata) / chunk_size)

    @staticmethod
    def _get_timestamp() -> str:
        return datetime.datetime.now().strftime("%H:%M:%S.%f")[:-3]

    def _print(self, str_to_print: str) -> None:
        print(f"* [{self._get_timestamp()}] {str_to_print}")

    def _validate_and_sanitize_input_data(
        self,
        adata: anndata.AnnData,
        feature_schema_name: str,
        count_matrix_name: str,
        feature_ids_column_name: str,
    ) -> anndata.AnnData:
        try:
            cas_feature_schema_list = self._feature_schemas_cache[feature_schema_name]
        except KeyError:
            cas_feature_schema_list = self.cas_api_service.get_cas_pca_002_schema_from_dump()
            self._feature_schemas_cache[feature_schema_name] = cas_feature_schema_list

        try:
            data_preparation.validate(
                adata=adata,
                cas_feature_schema_list=cas_feature_schema_list,
                feature_ids_column_name=feature_ids_column_name,
            )
        except exceptions.DataValidationError as e:
            if e.extra_features > 0:
                self._print(
                    f"The input data matrix has {e.extra_features} extra features compared to '{feature_schema_name}' "
                    f"CAS schema ({len(cas_feature_schema_list)}). "
                    f"Extra input features will be dropped."
                )
            if e.missing_features > 0:
                self._print(
                    f"The input data matrix has {e.missing_features} missing features compared to "
                    f"'{feature_schema_name}' CAS schema ({len(cas_feature_schema_list)}). "
                    f"Missing features will be imputed with zeros."
                )
            if e.extra_features == 0 and e.missing_features == 0:
                self._print(
                    f"Input datafile has all the necessary features as {feature_schema_name}, but it's still "
                    f"incompatible because of the different order. The features will be reordered "
                    f"according to {feature_schema_name}..."
                )
                self._print(
                    f"The input data matrix contains all of the features specified in '{feature_schema_name}' "
                    f"CAS schema but in a different order. The input features will be reordered according to "
                    f"'{feature_schema_name}'"
                )
            return data_preparation.sanitize(
                adata=adata,
                cas_feature_schema_list=cas_feature_schema_list,
                count_matrix_name=count_matrix_name,
                feature_ids_column_name=feature_ids_column_name,
            )
        else:
            self._print(f"The input data matrix conforms with the '{feature_schema_name}' CAS schema.")
            return adata

    async def _annotate_anndata_chunk(
        self,
        adata_bytes: bytes,
        results: t.List,
        chunk_index: int,
        chunk_start_i: int,
        chunk_end_i: int,
    ) -> None:
        """
        A wrapper around POST request that handles HTTP 500 and HTTP 401 status responses
        In case of HTTP 500 response resubmit task 2 times.
        In case of HTTP 501 print a message

        :param adata_bytes: Dumped bytes of `anndata.AnnData` chunk
        :param results: Results list that needs to be used to inplace the response from the server
        :param chunk_index: Consequent number of the chunk (e.g Chunk 1, Chunk 2)
        :param chunk_start_i: Index pointing to the main adata file start position of the current chunk
        :param chunk_end_i: Index pointing to the main adata file end position of the current chunk
        """
        number_of_cells = chunk_end_i - chunk_start_i
        for _ in range(self.num_attempts_per_chunk):
            try:
                results[chunk_index] = await self.cas_api_service.async_annotate_anndata_chunk(
                    adata_file_bytes=adata_bytes, number_of_cells=number_of_cells
                )

            except exceptions.HTTPError500:
                self._print(
                    f"Error occurred in CAS backend (HTTPError500), "
                    f"resubmitting chunk #{chunk_index + 1:2.0f} ({chunk_start_i:5.0f}, {chunk_end_i:5.0f}) to CAS ..."
                )
                pass
            except exceptions.HTTPError401:
                self._print("Unauthorized token. Please check your API token or request a new one.")
                break
            else:
                self._print(
                    f"Received the annotations for cell chunk #{chunk_index + 1:2.0f} ({chunk_start_i:5.0f}, {chunk_end_i:5.0f}) ..."
                )
                break

    async def _annotate_anndata_task(self, adata, results, chunk_size, start_time) -> None:
        """
        Submit chunks asynchronously as asyncio tasks

        :param adata: `anndata.AnnData` file to process
        :param results: Results list that is used to inplace the responses from the server
        :param chunk_size: Chunk size to split on
        :param start_time: Start time of the execution before this task was submitted
        """
        i, j = 0, chunk_size
        tasks = []
        number_of_chunks = self._get_number_of_chunks(adata, chunk_size=chunk_size)
        for chunk_index in range(number_of_chunks):
            chunk = adata[i:j, :]
            chunk_start_i = i
            chunk_end_i = i + len(chunk)
            self._print(
                f"Submitting cell chunk #{chunk_index + 1:2.0f} ({chunk_start_i:5.0f}, {chunk_end_i:5.0f}) to CAS ..."
            )

            chunk_bytes = _io.adata_to_bytes(adata=chunk)
            tasks.append(
                self._annotate_anndata_chunk(
                    adata_bytes=chunk_bytes,
                    results=results,
                    chunk_index=chunk_index,
                    chunk_start_i=chunk_start_i,
                    chunk_end_i=chunk_end_i,
                )
            )
            i = j
            j += chunk_size

        await asyncio.wait(tasks)

    def annotate_anndata(
        self,
        adata: "anndata.AnnData",
        chunk_size=2000,
        feature_schema_name: str = "default",
        count_matrix_name: str = "X",
        feature_ids_column_name: str = "index",
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Send an instance of :class:`anndata.AnnData` to the Cellarium Cloud backend for annotations. The function
        splits the ``adata`` into smaller chunks and asynchronously sends them to the backend API service. Each chunk is
        of equal size, except for the last one, which may be smaller. The backend processes these chunks in parallel.

        :param adata: :class:`anndata.AnnData` instance to annotate
        :param chunk_size: Size of chunks to split on
        :param feature_schema_name: feature schema name to use for data preparation. |br|
            `Allowed Values:` Feature schema name from the :attr:`feature_schemas` list or ``"default"``
            keyword, which refers to the default selected feature schema in the Cellarium backend. |br|
            `Default:` ``"default"``
        :param count_matrix_name:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choice of either ``"X"``  or ``"raw.X"`` in order to use ``adata.X`` or ``adata.raw.X``,
             respectively |br|
            `Default:` ``"X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :return: A list of dictionaries with annotations for each of the cells from input adata
        """
        start = time.time()
        self._print("CAS v1.0 (Model ID: PCA_002)")
        self._print(f"Total number of input cells: {len(adata)}")
        assert feature_schema_name == "default" or feature_schema_name in self.feature_schemas, (
            "`feature_schema_name` should have a value of either 'default' or one of the values from "
            "`feature_schemas`."
        )
        adata = self._validate_and_sanitize_input_data(
            adata=adata,
            feature_schema_name=feature_schema_name,
            count_matrix_name=count_matrix_name,
            feature_ids_column_name=feature_ids_column_name,
        )
        # TODO: Deprecation of having `raw` is needed in newer versions of CAS.
        adata.raw = adata
        number_of_chunks = math.ceil(len(adata) / chunk_size)
        results = [[] for _ in range(number_of_chunks)]
        loop = asyncio.get_event_loop()
        task = loop.create_task(
            self._annotate_anndata_task(adata=adata, results=results, chunk_size=chunk_size, start_time=start)
        )
        loop.run_until_complete(task)
        result = functools.reduce(operator.iconcat, results, [])
        self._print(f"Total wall clock time: {f'{time.time() - start:10.4f}'} seconds")
        self._print("Finished!")
        return result

    def annotate_anndata_file(
        self,
        filepath: str,
        chunk_size=2000,
        feature_schema_name: str = "default",
        count_matrix_name: str = "X",
        feature_ids_column_name: str = "index",
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Read the 'h5ad' file into a :class:`anndata.AnnData` matrix and apply the :meth:`annotate_anndata` method to it.

        :param filepath: Filepath of the local :class:`anndata.AnnData` matrix
        :param chunk_size: Size of chunks to split on
        :param feature_schema_name: feature schema name to use for data preparation. |br|
            `Allowed Values:` Feature schema name from the :attr:`feature_schemas` list or ``"default"``
            keyword, which refers to the default selected feature schema in the Cellarium backend. |br|
            `Default:` ``"default"``
        :param count_matrix_name:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choice of either ``"X"``  or ``"raw.X"`` in order to use ``adata.X`` or ``adata.raw.X``,
             respectively |br|
            `Default:` ``"X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :return: A list of dictionaries with annotations for each of the cells from input adata
        """
        adata = anndata.read_h5ad(filename=filepath)
        return self.annotate_anndata(
            adata=adata,
            chunk_size=chunk_size,
            feature_schema_name=feature_schema_name,
            count_matrix_name=count_matrix_name,
            feature_ids_column_name=feature_ids_column_name,
        )

    def annotate_10x_h5_file(
        self,
        filepath: str,
        chunk_size: int = 2000,
        feature_schema_name: str = "default",
        count_matrix_name: str = "X",
        feature_ids_column_name: str = "index",
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Parse the 10x 'h5' matrix and apply the :meth:`annotate_anndata` method to it.

        :param filepath: Filepath of the local 'h5' matrix
        :param chunk_size: Size of chunks to split on
        :param feature_schema_name: feature schema name to use for data preparation. |br|
            `Allowed Values:` Feature schema name from the :attr:`feature_schemas` list or ``"default"``
            keyword, which refers to the default selected feature schema in the Cellarium backend. |br|
            `Default:` ``"default"``
        :param count_matrix_name:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choice of either ``"X"``  or ``"raw.X"`` in order to use ``adata.X`` or ``adata.raw.X``,
             respectively |br|
            `Default:` ``"X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :return: A list of dictionaries with annotations for each of the cells from input adata
        """
        adata = _io.read_10x_h5(filepath)

        return self.annotate_anndata(
            adata=adata,
            chunk_size=chunk_size,
            feature_schema_name=feature_schema_name,
            count_matrix_name=count_matrix_name,
            feature_ids_column_name=feature_ids_column_name,
        )
