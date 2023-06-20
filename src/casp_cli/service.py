import asyncio
import datetime
import functools
import math
import operator
import os
import time
import typing as t

import aiohttp
import nest_asyncio

from casp_cli import _read_data, exceptions

if t.TYPE_CHECKING:
    import anndata

nest_asyncio.apply()


class _BaseService:
    BACKEND_URL: str

    def __init__(self, api_token: str, *args, **kwargs):
        """
        Base clas for communicating with a Cellarium Cloud API service
        It leverages async request library `aiohttp` to asynchrounously execute HTTP requests.

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

    async def post(
        self, endpoint: str, file, data: t.Optional[t.Dict] = None, headers: t.Optional[t.Dict] = None
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
                if resp.status == 401:
                    raise exceptions.HTTPError401
                elif resp.status == 403:
                    raise exceptions.HTTPError403
                elif resp.status == 500:
                    raise exceptions.HTTPError500

                return await resp.json()


class CASClientService(_BaseService):
    """
    Service that is designed to call Cellarium Cloud Backend API service
    """

    BACKEND_URL = "https://cas-manager-vi7nxpvk7a-uc.a.run.app"

    def _get_number_of_chunks(self, adata, chunk_size):
        return math.ceil(len(adata) / chunk_size)

    def _get_timestamp(self) -> str:
        return datetime.datetime.now().strftime("%H:%M:%S.%f")[:-3]

    def _print(self, str_to_print: str) -> None:
        print(f"* [{self._get_timestamp()}] {str_to_print}")

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
                results[chunk_index] = await self.post("annotate", file=adata_bytes, data=data)
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
        :param chunk_sie: Chunk size to split on
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

    def annotate_anndata(self, adata: "anndata.AnnData", chunk_size=2000) -> t.List[t.Dict[str, t.Any]]:
        """
        Send an anndata object to Cellarium Cloud backend to get annotations. Split the input
        adata to smaller chunks and send them all asynchronously to the backend API service.
        All of the chunks have equal size apart of the last one, which is usually smaller than
        the rest of the chunks. Backend API processes all of these chunks in parallel and
        returns them as soon as they are ready.

        :param adata: AnnData instance to annotate
        :param chunk_size: Size of chunks to split on
        :return: A list of dictionaries with annotations for each of the cells from input adata
        """
        start = time.time()
        self._print("CAS v1.0 (Model ID: PCA_002)")
        self._print(f"Total number of input cells: {len(adata)}")
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

    def annotate_10x_h5(self, filepath: str, chunk_size=2000) -> t.List[t.Dict[str, t.Any]]:
        adata = _read_data.read_10x_h5(filepath)
        return self.annotate_anndata(adata=adata, chunk_size=chunk_size)
