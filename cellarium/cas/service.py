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
    def raise_response_exception(status_code: int, detail: str) -> None:
        """
        Raise an exception based on the status code returned by the server, including the detail message

        :param status_code: HTTP status code
        :param detail: Detail message returned by the server
        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError
        """
        message = f"Server returned status code {status_code}, Detail: {detail}"
        if status_code == 401:
            raise exceptions.HTTPError401(message)
        elif status_code == 403:
            raise exceptions.HTTPError403(message)
        elif status_code == 500:
            raise exceptions.HTTPError500(message)
        else:
            raise exceptions.HTTPError(message)

    def get(self, endpoint: str) -> requests.Response:
        """
        Make a GET request to backend service

        :param endpoint: Endpoint string without a leading slash

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError

        :return: Response object
        """
        url = self._get_endpoint_url(endpoint)
        headers = {"Authorization": f"Bearer {self.api_token}"}
        response = requests.get(url=url, headers=headers)

        status_code = response.status_code
        if status_code < 200 or status_code >= 300:
            try:
                response_detail = response.json()["detail"]
            except (json.decoder.JSONDecodeError, KeyError):
                response_detail = response.text

            self.raise_response_exception(status_code=status_code, detail=response_detail)

        return response

    def get_json(self, endpoint: str) -> t.Union[t.Dict, t.List]:
        """
        Make a GET request to backend service and return JSON response

        :param endpoint: Endpoint string without a leading slash

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError

        :return: JSON response
        """
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

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError

        :return: JSON response
        """
        url = self._get_endpoint_url(endpoint=endpoint)
        _data = {}
        _headers = {"Authorization": f"Bearer {self.api_token}"}

        if data is not None:
            _data.update(**data)
        if headers is not None:
            _headers.update(**headers)

        form_data = aiohttp.FormData()
        form_data.add_field("file", file, filename="adata.h5ad")

        for key, value in data.items():
            form_data.add_field(key, value)

        ssl_context = ssl.create_default_context(cafile=certifi.where())
        conn = aiohttp.TCPConnector(ssl=ssl_context)

        async with aiohttp.ClientSession(connector=conn) as session:
            async with session.post(url, data=form_data, headers=_headers) as response:
                status_code = response.status
                if status_code < 200 or status_code >= 300:
                    try:
                        response_body = await response.json()
                        response_detail = response_body["detail"]
                    except (json.decoder.JSONDecodeError, KeyError):
                        response_detail = await response.text()

                    self.raise_response_exception(status_code=status_code, detail=response_detail)

                return await response.json()


class CASAPIService(_BaseService):
    """
    Class with all the API methods of Cellarium Cloud CAS infrastructure.
    """

    BACKEND_URL = "https://cas-api-1-3-xdev-vi7nxpvk7a-uc.a.run.app"

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

    def get_model_list(self) -> t.List[t.Dict[str, t.Any]]:
        """
        Retrieve list of all models that are in CAS

        Refer to API Docs: {BACKEND_URL}/docs#/default/list_models_list_models_get

        :return: List of models
        """
        return self.get_json(endpoint=endpoints.LIST_MODELS)

    async def async_annotate_anndata_chunk(
        self, adata_file_bytes: t.ByteString, number_of_cells: int, model_name: str, include_dev_metadata: bool
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Request Cellarium Cloud infrastructure to annotate an input anndata file

        Refer tp API Docs: {BACKEND_URL}/docs#/default/annotate_annotate_post

        :param adata_file_bytes: Validated anndata file
        :param number_of_cells: Number of cells being processed in this dataset
        :param model_name: Name of the model to use.
        :param include_dev_metadata: Whether to include dev metadata in the response.
        :return: A list of dictionaries with annotations.
        """
        request_data = {
            "number_of_cells": str(number_of_cells),
            "model_name": model_name,
            "include_dev_metadata": str(include_dev_metadata),
        }
        return await self.async_post(endpoints.ANNOTATE, file=adata_file_bytes, data=request_data)
