import json
import ssl
import typing as t

import aiohttp
import certifi
import nest_asyncio
import requests

from cellarium.cas import endpoints, exceptions

nest_asyncio.apply()


class _BaseService:
    BACKEND_URL: str

    def __init__(self, api_token: str, *args, **kwargs):
        """
        Base class for communicating with a Cellarium Cloud API service
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

    def get_json(self, endpoint: str) -> t.Union[t.Dict, t.List]:
        return self.get(endpoint=endpoint).json()

    async def async_post(
        self,
        endpoint: str,
        file,
        data: t.Optional[t.Dict] = None,
        headers: t.Optional[t.Dict] = None,
    ) -> t.Union[t.List, t.Dict]:
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

        ssl_context = ssl.create_default_context(cafile=certifi.where())
        conn = aiohttp.TCPConnector(ssl=ssl_context)

        async with aiohttp.ClientSession(connector=conn) as session:
            async with session.post(url, data=form_data, headers=_headers) as resp:
                self.__validate_response_code(resp.status)
                return await resp.json()


class CASAPIService(_BaseService):
    """
    Class with all the API methods of Cellarium Cloud CAS infrastructure.
    """

    BACKEND_URL = "https://cas-api-test-2-vi7nxpvk7a-uc.a.run.app"

    def validate_token(self) -> None:
        """
        Validate user given API token.
        Would raise 401 Unauthorized if token is invalid.

        Refer to API Docs: {BACKEND_URL}/docs#/default/validate_token_validate_token_get

        :return: Void
        """
        self.get(endpoint=endpoints.VALIDATE_TOKEN)

    def get_application_info(self) -> t.Dict[str, str]:
        """
        Retrieve General Application Info. This Includes default schema, version, model information, etc.

        Refer to API Docs: {BACKEND_URL}/docs#/default/application_info_application_info_get

        :return: Dictionary with application info
        """
        return self.get_json(endpoint=endpoints.APPLICATION_INFO)

    def get_feature_schemas(self) -> t.List[str]:
        """
        Retrieve a list of feature schemas that exist in Cellarium Cloud CAS

        Refer to API Docs: {BACKEND_URL}/docs#/default/get_feature_schemas_feature_schemas_get

        :return: List of feature schema names
        """
        return [x["schema_name"] for x in self.get_json(endpoint=endpoints.GET_FEATURE_SCHEMAS)]

    def get_feature_schema_by(self, name: str) -> t.List[str]:
        """
        Retrieve feature schema by name

        Refer to API Docs: {BACKEND_URL}/docs#/default/get_feature_schema_by_feature_schema__schema_name__get

        :param name: Name of feature schema
        :return: List of feature ids
        """
        return self.get_json(endpoint=endpoints.GET_SCHEMA_BY_NAME.format(schema_name=name))

    @staticmethod
    def get_cas_pca_002_schema_from_dump() -> t.List[str]:
        """
        This method should be deprecated and used only before Cellarium Cloud CAS backend will have feature schemas
        API methods
        :return: list with gene ids
        """
        import os

        curr_path = os.path.dirname(os.path.realpath(__file__))

        with open(f"{curr_path}/assets/cellarium_cas_tx_pca_002_grch38_2020_a.json", "r") as f:
            cas_feature_schema_list = json.loads(f.read())

        return cas_feature_schema_list

    async def async_annotate_anndata_chunk(
        self, adata_file_bytes: t.ByteString, number_of_cells: int
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Request Cellarium Cloud infrastructure to annotate an input anndata file

        Refer tp API Docs: {BACKEND_URL}/docs#/default/annotate_annotate_post

        :param adata_file_bytes: Validated anndata file
        :param number_of_cells: Number of cells being processed in this dataset
        :return: A list of dictionaries with annotations.
        """
        request_data = {"number_of_cells": str(number_of_cells)}
        return await self.async_post(endpoints.ANNOTATE, file=adata_file_bytes, data=request_data)
