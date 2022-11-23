import asyncio
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

    def __init__(self):
        self.session = aiohttp.ClientSession(raise_for_status=True)

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

    async def _annotate_anndata_chunk(self, adata_bytes, results, chunk_index, start_time) -> None:
        data = {"myfile": adata_bytes}
        results[chunk_index] = await self.post("annotate", data=data)
        print(f"--- FINISHED PROCESSING CHUNK {time.time() - start_time} ---")

    async def _annotate_anndata_task(self, adata, results, chunk_size, start_time) -> None:
        i, j = 0, chunk_size
        tasks = []
        number_of_chunks = self._get_number_of_chunks(adata, chunk_size=chunk_size)
        for chunk_index in range(number_of_chunks):
            chunk = adata[i:j, :]
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
                    start_time=start_time,
                )
            )
            i = j
            j += chunk_size

        await asyncio.wait(tasks)

    def annotate_anndata(self, adata: "anndata.AnnData", chunk_size=2000) -> t.List:
        number_of_chunks = math.ceil(len(adata) / chunk_size)
        results = [None] * number_of_chunks
        start = time.time()
        loop = asyncio.get_event_loop()
        task = loop.create_task(
            self._annotate_anndata_task(adata=adata, results=results, chunk_size=chunk_size, start_time=start)
        )
        loop.run_until_complete(task)
        print(f"--- TOTAL TIME SPENT {time.time() - start} -- ")
        return functools.reduce(operator.iconcat, results, [])
