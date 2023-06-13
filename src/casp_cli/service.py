import asyncio
import datetime
import functools
import math
import operator
import os
import time
import typing as t

import aiohttp
import anndata
import nest_asyncio
import requests

from casp_cli import data_preparation, endpoints, exceptions

nest_asyncio.apply()


class _BaseService:
    BACKEND_URL: str

    def __init__(self, api_token: str, *args, **kwargs):
        """
        Base clas for communicating with a Cellarium Cloud API service
        It leverages async request library `aiohttp` to asynchronously execute HTTP requests.

        :param api_token: A token that could be authenticated by Cellarium Cloud Backend API service
        """
        self.api_token = api_token
        super().__init__(*args, **kwargs)

    @classmethod
    def _get_endpoint_url(cls, endpoint: str) -> str:
        """
        Configure a specific method endpoint from backend url and endpoint

        :param endpoint: Endpoint string without a leading slash
        :return: Full url with backend domains/subdomains and endpoint joint
        """
        return f"{cls.BACKEND_URL}/{endpoint}"

    @staticmethod
    def __validate_response_code(response_code):
        if response_code == 401:
            raise exceptions.HTTPError401
        elif response_code == 403:
            raise exceptions.HTTPError403
        elif response_code == 500:
            raise exceptions.HTTPError500

    def get(self, endpoint: str) -> requests.Response:
        url = self._get_endpoint_url(endpoint)
        headers = {"Authorization": f"Bearer {self.api_token}"}
        response = requests.get(url=url, headers=headers)
        self.__validate_response_code(response.status_code)
        return response

    def get_json(self, endpoint: str) -> dict:
        return self.get(endpoint=endpoint).json()

    async def async_post(
        self,
        endpoint: str,
        file,
        data: t.Optional[t.Dict] = None,
        headers: t.Optional[t.Dict] = None,
    ) -> t.Dict:
        """
        Make an async multiform POST request to backend service

        :param endpoint: Endpoint string without a leading slash
        :param file: Byte file to attach to request
        :param data: Dictionary to include to HTTP POST request body
        :param headers: Dictionary to include to HTTP POST request Headers
        """
        url = self._get_endpoint_url(endpoint=endpoint)
        _data = {}
        _headers = {"Authorization": f"Bearer {self.api_token}"}

        if data is not None:
            _data.update(**data)
        if headers is not None:
            _headers.update(**headers)

        form_data = aiohttp.FormData()
        form_data.add_field("myfile", file, filename="adata.h5ad")

        for key, value in data.items():
            form_data.add_field(key, value)

        async with aiohttp.ClientSession() as session:
            async with session.post(url, data=form_data, headers=_headers) as resp:
                self.__validate_response_code(resp.status)
                return await resp.json()


class CASClientService(_BaseService):
    """
    Service that is designed to call Cellarium Cloud Backend API service
    """

    BACKEND_URL = "https://cas-api-fg-gene-schema-2-vi7nxpvk7a-uc.a.run.app"

    def __init__(self, api_token: str, *args, **kwargs):
        super().__init__(api_token=api_token, *args, **kwargs)
        self._print("Connecting with Cellarium Cloud backend...")
        # Validate Auth Token
        self.get(endpoint=endpoints.VALIDATE_TOKEN)
        # Get application info
        app_info = self.get_json(endpoint=endpoints.APPLICATION_INFO)
        self.default_feature_schema_name = app_info["default_feature_schema"]
        # Retrieve feature schemas
        self.feature_schemas = [x["schema_name"] for x in self.get_json(endpoint=endpoints.GET_FEATURE_SCHEMAS)]
        self._app_version = app_info["application_version"]
        self._feature_schemas_cache = {}
        self._print(
            (
                f"Authorized. CAS v {self._app_version}. You're using a default feature schema "
                f"{self.default_feature_schema_name}."
            )
        )

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
            cas_feature_schema_list = self.get_json(
                endpoint=endpoints.GET_SCHEMA_BY_NAME.format(schema_name=feature_schema_name)
            )
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
                    f"Input datafile has {e.extra_features} more features than {feature_schema_name} schema in CAS. "
                    f"They will be omitted..."
                )
            if e.missing_features > 0:
                self._print(
                    f"Input datafile has {e.missing_features} missing features compared with {feature_schema_name} in CAS. "
                    f"They will be filled with zeros..."
                )
            if e.extra_features == 0 and e.missing_features == 0:
                self._print(
                    f"Input datafile has all the necessary features as {feature_schema_name}, but it's still "
                    f"incompatible because of the different order. The features will be reordered "
                    f"according to {feature_schema_name}..."
                )
            return data_preparation.sanitize(
                adata=adata,
                cas_feature_schema_list=cas_feature_schema_list,
                count_matrix_name=count_matrix_name,
                feature_ids_column_name=feature_ids_column_name,
            )
        else:
            self._print(f"Input dataset is validated and fully compatible with {feature_schema_name} schema...")
            return adata

    async def _annotate_anndata_chunk(
        self,
        adata_bytes: bytes,
        results: t.List,
        chunk_index: int,
        chunk_start_i: int,
        chunk_end_i: int,
        start_time: float,
    ) -> None:
        """
        A wrapper around POST request that handles HTTP 500 and HTTP 401 status responses
        In case of HTTP 500 response resubmit task 2 times.
        In case of HTTP 501 print a message

        :param adata_bytes: Dumped bytes of anndata chunk
        :param results: Results list that needs to be used to inplace the response from the server
        :param chunk_index: Consequent number of the chunk (e.g Chunk 1, Chunk 2)
        :param chunk_start_i: Index pointing to the main adata file start position of the current chunk
        :param chunk_end_i: Index pointing to the main adata file end position of the current chunk
        :param start_time: Start time of the execution before any task was submitted
        """
        number_of_cells = chunk_end_i - chunk_start_i
        data = {"number_of_cells": str(number_of_cells)}
        for _ in range(3):
            try:
                results[chunk_index] = await self.async_post(endpoints.ANNOTATE, file=adata_bytes, data=data)
            except exceptions.HTTPError500:
                self._print(
                    f"Something went wrong, resubmitting chunk #{chunk_index + 1:2.0f} ({chunk_start_i:5.0f}, {chunk_end_i:5.0f}) to CAS ..."
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

        :param adata: AnnData file to process
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
            tmp_file_name = f"chunk_{chunk_index}.h5ad"
            chunk.write(tmp_file_name, compression="gzip")
            with open(tmp_file_name, "rb") as f:
                chunk_bytes = f.read()

            os.remove(tmp_file_name)
            tasks.append(
                self._annotate_anndata_chunk(
                    adata_bytes=chunk_bytes,
                    results=results,
                    chunk_index=chunk_index,
                    chunk_start_i=chunk_start_i,
                    chunk_end_i=chunk_end_i,
                    start_time=start_time,
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
        Send an anndata object to Cellarium Cloud backend to get annotations. Split the input
        adata to smaller chunks and send them all asynchronously to the backend API service.
        All the chunks have equal size apart of the last one, which is usually smaller than
        the rest of the chunks. Backend API processes all of these chunks in parallel and
        returns them as soon as they are ready.

        :param adata: AnnData instance to annotate
        :param chunk_size: Size of chunks to split on
        :param feature_schema_name: feature schema name to use for data preparation. if 'default' default schema will
            be used.
        :param count_matrix_name:  Where to obtain a feature expression count matrix from. Choice of: 'X', 'raw.X'
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. Default `index`.
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
