import asyncio
import json
import ssl
import typing as t

import aiohttp
import certifi
import nest_asyncio
import requests
from aiohttp import client_exceptions

from cellarium.cas import endpoints, exceptions

nest_asyncio.apply()

AIOHTTP_TOTAL_TIMEOUT_SECONDS = 220
AIOHTTP_READ_TIMEOUT_SECONDS = 200


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
        elif status_code == 500 or status_code == 502 or status_code == 503 or status_code == 504:
            raise exceptions.HTTPError5XX(message)
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

    async def _aiohttp_async_post(
        self, url: str, form_data: t.Optional[aiohttp.FormData] = None, headers: t.Optional[t.Dict[str, t.Any]] = None
    ):
        """
        Make an async POST request to backend service with timeout and SSL context, handle exceptions and return JSON.

        Create custom SSL context with `certifi` package to avoid SSL errors. Use :class:`aiohttp.ClientTimeout` to set
        timeout for the request. Handle exceptions and raise custom exceptions.

        :param url: Endpoint URL
        :param form_data: :class:`aiohttp.FormData` object
        :param headers: Headers to include in the request

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError, HTTPClientError

        :return: JSON response
        """
        timeout = aiohttp.ClientTimeout(total=AIOHTTP_TOTAL_TIMEOUT_SECONDS, sock_read=AIOHTTP_READ_TIMEOUT_SECONDS)
        ssl_context = ssl.create_default_context(cafile=certifi.where())
        connector = aiohttp.TCPConnector(ssl=ssl_context)

        async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
            try:
                async with session.post(url, data=form_data, headers=headers) as response:
                    status_code = response.status
                    if status_code < 200 or status_code >= 300:
                        try:
                            response_body = await response.json()
                            response_detail = response_body["detail"]
                        except (json.decoder.JSONDecodeError, client_exceptions.ClientResponseError):
                            response_detail = await response.text()
                        except KeyError:
                            print("Response body doesn't have a 'detail' key, returning full response body")
                            response_detail = str(await response.json())

                        self.raise_response_exception(status_code=status_code, detail=response_detail)
                    try:
                        return await response.json()
                    except client_exceptions.ClientPayloadError as e:
                        raise exceptions.HTTPClientError(f"Failed to parse response body: {e.__class__.__name__}: {e}")

            except (client_exceptions.ServerTimeoutError, asyncio.TimeoutError) as e:
                raise exceptions.HTTPClientError(f"Client Timeout Error: {e}")
            except (client_exceptions.ClientConnectionError, client_exceptions.ClientResponseError) as e:
                raise exceptions.HTTPClientError(f"Client Connection Error: {e.__class__.__name__}: {e}")
            except client_exceptions.ClientError as e:
                raise exceptions.HTTPClientError(f"Unknown Error: {e.__class__.__name__}: {e}")

    async def async_post(
        self,
        endpoint: str,
        file: t.ByteString,
        data: t.Optional[t.Dict] = None,
        headers: t.Optional[t.Dict] = None,
    ) -> t.Union[t.List, t.Dict]:
        """
        Make an async multiform POST request to backend service. Include adata file as a byte file in the request body,
        and authorization token in the headers.

        For details about custom settings for SSL context and timeout, refer to :meth:`_aiohttp_async_post` method.

        :param endpoint: Endpoint string without a leading slash
        :param file: Byte file to attach to request
        :param data: Dictionary to include to HTTP POST request body
        :param headers: Dictionary to include to HTTP POST request Headers

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError, HTTPClientError

        :return: JSON response
        """
        url = self._get_endpoint_url(endpoint=endpoint)
        _data = data if data is not None else {}
        _headers = {"Authorization": f"Bearer {self.api_token}"}

        if headers is not None:
            _headers.update(headers)

        form_data = aiohttp.FormData()
        form_data.add_field("file", file, filename="adata.h5ad")

        for key, value in _data.items():
            form_data.add_field(key, value)

        return await self._aiohttp_async_post(url=url, form_data=form_data, headers=_headers)


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
