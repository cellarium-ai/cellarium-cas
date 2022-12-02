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

if t.TYPE_CHECKING:
    import anndata
nest_asyncio.apply()


class _BaseService:
    BACKEND_URL: str

    @classmethod
    def _get_endpoint_url(cls, endpoint: str) -> str:
        return f"{cls.BACKEND_URL}/{endpoint}"

    async def post(self, endpoint: str, data) -> t.Dict:
        url = self._get_endpoint_url(endpoint=endpoint)

        async with aiohttp.ClientSession() as session:
            async with session.post(url, data=data) as resp:
                return await resp.json()


class CASPClientService(_BaseService):
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
        data = {"myfile": adata_bytes}
        for _ in range(3):
            try:
                results[chunk_index] = await self.post("annotate", data=data)
            except Exception:
                self._print(f"Something went wrong, resubmitting chunk #{chunk_index + 1:2.0f} ({chunk_start_i:5.0f}, {chunk_end_i:5.0f}) to CAS ...")
                pass
            else:
                break
        self._print(f"Received the annotations for cell chunk #{chunk_index + 1:2.0f} ({chunk_start_i:5.0f}, {chunk_end_i:5.0f}) ...")

    async def _annotate_anndata_task(self, adata, results, chunk_size, start_time) -> None:
        i, j = 0, chunk_size
        tasks = []
        number_of_chunks = self._get_number_of_chunks(adata, chunk_size=chunk_size)
        for chunk_index in range(number_of_chunks):
            chunk = adata[i:j, :]
            chunk_start_i = i
            chunk_end_i = i + len(chunk)
            self._print(f"Submitting cell chunk #{chunk_index + 1:2.0f} ({chunk_start_i:5.0f}, {chunk_end_i:5.0f}) to CAS ...")
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

    def annotate_anndata(self, adata: "anndata.AnnData", chunk_size=2000) -> t.List:
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
