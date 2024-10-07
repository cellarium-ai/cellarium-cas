import asyncio
import json
import ssl
import typing as t
from contextlib import ContextDecorator
from contextvars import ContextVar
from uuid import UUID, uuid4

import aiohttp
import certifi
import nest_asyncio
import requests
from aiohttp import client_exceptions

from cellarium.cas import constants, endpoints, exceptions, settings
from cellarium.cas.logging import logger

if settings.is_interactive_environment():
    logger.debug("Running in an interactive environment, applying nest_asyncio")
    nest_asyncio.apply()

# Context variable to track the action id for the current context
client_action_id = ContextVar("action_id", default=None)


class action_context_manager(ContextDecorator):
    """
    Handle the lifecycle of an action id for the current context.

    Add the `@action_context_manager` decorator around methods in the client methods that should
    have a common action id.
    """

    def __enter__(self):
        self.t = client_action_id.set(uuid4())
        return self

    def __exit__(self, *exc):
        client_action_id.reset(self.t)


class _BaseService:
    """
    Base class for communicating with a Cellarium Cloud API service
    It leverages async request library `aiohttp` to asynchronously execute HTTP requests.

    :param api_token: A token that could be authenticated by Cellarium Cloud Backend API service
    :param api_url: URL of the Cellarium Cloud Backend API service
    :param client_session_id: A unique identifier for the current client session
    """

    def __init__(
        self,
        api_token: str,
        api_url: str = settings.CELLARIUM_CLOUD_BACKEND_URL,
        client_session_id: UUID = None,
        *args,
        **kwargs,
    ):
        self.api_token = api_token
        self.api_url = api_url
        self.client_session_id = client_session_id
        super().__init__(*args, **kwargs)

    def _get_endpoint_url(self, endpoint: str) -> str:
        """
        Configure a specific method endpoint from backend url and endpoint

        :param endpoint: Endpoint string without a leading slash
        :return: Full url with backend domains/subdomains and endpoint joined
        """
        return f"{self.api_url}/{endpoint}"

    def _get_headers(self) -> t.Dict[str, str]:
        """
        Get the headers to include in the request

        :return: Headers dictionary
        """
        headers = {constants.Headers.authorization: f"Bearer {self.api_token}"}
        if self.client_session_id:
            headers[constants.Headers.client_session_id] = str(self.client_session_id)
        action_id = client_action_id.get()
        if action_id:
            headers[constants.Headers.client_action_id] = str(action_id)
        return headers

    @staticmethod
    def raise_response_exception(status_code: int, detail: str) -> None:
        """
        Raise an exception based on the status code returned by the server, including the detail message

        :param status_code: HTTP status code
        :param detail: Detail message returned by the server
        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError
        """
        message = f"Server returned status code {status_code}, Detail: {detail}"
        if status_code == constants.HTTP.STATUS_401_UNAUTHORIZED:
            raise exceptions.HTTPError401(message)
        elif status_code == constants.HTTP.STATUS_403_FORBIDDEN:
            raise exceptions.HTTPError403(message)
        elif status_code == constants.HTTP.STATUS_NOT_FOUND:
            raise exceptions.HTTPError404(message)

        elif (
            constants.HTTP.STATUS_500_INTERNAL_SERVER_ERROR
            <= status_code
            <= constants.HTTP.STATUS_511_NETWORK_AUTHENTICATION_REQUIRED
        ):
            raise exceptions.HTTPError5XX(message)
        else:
            raise exceptions.HTTPError(message)

    def __validate_requests_response(self, response: requests.Response) -> None:
        """
        Validate requests response and raise an exception if response status code is not 200

        :param response: Response object

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError
        """
        status_code = response.status_code
        if not (constants.HTTP.STATUS_200_OK <= status_code <= constants.HTTP.STATUS_226_IM_USED):
            # When response status code is not 2XX
            try:
                response_detail = response.json()["detail"]
            except (json.decoder.JSONDecodeError, KeyError):
                response_detail = response.text

            self.raise_response_exception(status_code=status_code, detail=response_detail)

    def get(self, endpoint: str) -> requests.Response:
        """
        Make a GET request to backend service

        :param endpoint: Endpoint string without a leading slash

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError

        :return: Response object
        """
        url = self._get_endpoint_url(endpoint)
        headers = self._get_headers()
        response = requests.get(url=url, headers=headers)
        self.__validate_requests_response(response=response)
        return response

    def post(self, endpoint: str, data: t.Optional[t.Dict] = None) -> requests.Response:
        url = self._get_endpoint_url(endpoint)
        headers = self._get_headers()
        response = requests.post(url=url, headers=headers, json=data)
        self.__validate_requests_response(response=response)
        return response

    def get_json(self, endpoint: str) -> t.Union[t.Dict, t.List]:
        """
        Make a GET request to backend service and return JSON response

        :param endpoint: Endpoint string without a leading slash

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError

        :return: JSON response
        """
        return self.get(endpoint=endpoint).json()

    def post_json(self, endpoint: str, data: t.Optional[t.Dict] = None) -> t.Union[t.Dict, t.List]:
        """
        Make a POST request to backend service and return JSON response

        :param endpoint: Endpoint string without a leading slash
        :param data: Dictionary to include to HTTP POST request body

        :raises: HTTPError401, HTTPError403, HTTPError500, HTTPBaseError

        :return: JSON response
        """
        return self.post(endpoint=endpoint, data=data).json()

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
        timeout = aiohttp.ClientTimeout(
            total=settings.AIOHTTP_TOTAL_TIMEOUT_SECONDS, sock_read=settings.AIOHTTP_READ_TIMEOUT_SECONDS
        )
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
                            logger.warning("Response body doesn't have a 'detail' key, returning full response body")
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
        _headers = self._get_headers()

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

    def validate_token(self) -> t.Dict[str, str]:
        """
        Validate user given API token.
        Would raise 401 Unauthorized if token is invalid.

        Refer to API Docs:
        {api_url}/api/docs#/cellarium-general/validate_token_api_cellarium_general_validate_token_get

        :return: Dictionary with basic user information
        """
        return self.get_json(endpoint=endpoints.VALIDATE_TOKEN)

    def get_application_info(self) -> t.Dict[str, str]:
        """
        Retrieve General Application Info. This Includes default schema, version, model information, etc.

        Refer to API Docs:
        {api_url}/api/docs#/cellarium-general/application_info_api_cellarium_general_application_info_get

        :return: Dictionary with application info
        """
        return self.get_json(endpoint=endpoints.APPLICATION_INFO)

    def feedback_opt_out(self) -> t.Dict[str, str]:
        """
        Opt out of future feedback requests.

        Refer to API Docs:
        {api_url}/api/docs#/cellarium-general/feedback/opt-out
        :return: Dictionary with basic user information
        """
        return self.post_json(endpoint=endpoints.FEEDBACK_OPT_OUT)

    def get_feedback_answer_link(self) -> str:
        """
        Retrieve a link to answer feedback questions

        Refer to API Docs:
        {api_url}/api/docs#/cellarium-general/feedback/answer
        :return: Link to answer feedback questions
        """
        return self._get_endpoint_url(
            endpoints.FEEDBACK_ANSWER.format(
                client_session_id=self.client_session_id, client_action_id=client_action_id.get()
            )
        )

    def validate_version(self, version_str: str) -> t.Dict[str, t.Any]:
        """
        Validate client version with the server to see if it is compatible.
        Would raise 400 Bad Request if version is not compatible.

        :return: Void
        """
        return self.post_json(endpoint=endpoints.VALIDATE_VERSION, data={"client_version": version_str})

    def get_feature_schemas(self) -> t.List[str]:
        """
        Retrieve a list of feature schemas that exist in Cellarium Cloud CAS

        Refer to API Docs:
        {api_url}/api/docs#/cellarium-general/get_feature_schemas_api_cellarium_general_feature_schemas_get

        :return: List of feature schema names
        """
        return [x["schema_name"] for x in self.get_json(endpoint=endpoints.GET_FEATURE_SCHEMAS)]

    def get_feature_schema_by(self, name: str) -> t.List[str]:
        """
        Retrieve feature schema by name

        Refer to API Docs:
        {api_url}/api/docs#/cellarium-general/get_feature_schema_by_api_cellarium_general_feature_schema__schema_name__get

        :param name: Name of feature schema

        :return: List of feature ids
        """
        return self.get_json(endpoint=endpoints.GET_SCHEMA_BY_NAME.format(schema_name=name))

    def get_model_list(self) -> t.List[t.Dict[str, t.Any]]:
        """
        Retrieve list of all models that are in CAS

        Refer to API Docs:
        {api_url}/api/docs#/cellarium-general/get_model_list_api_cellarium_general_list_models_get

        :return: List of models
        """
        return self.get_json(endpoint=endpoints.LIST_MODELS)

    def get_user_quota(self) -> t.Dict[str, t.Any]:
        """
        Retrieve user quota information

        Refer to API Docs:
        {BACKEND_URL}/api/docs

        :return: User quota information
        """
        return self.get_json(endpoint=endpoints.GET_USER_QUOTA)

    def query_cells_by_ids(
        self, cell_ids: t.List[int], metadata_feature_names: t.List[str]
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Retrieve cells by their ids from Cellarium Cloud database.

        Refer to API Docs:
        {api_url}/api/docs#/cell-analysis/get_cells_by_ids_api_cellarium_cas_query_cells_by_ids_post

        :param cell_ids: List of cell ids from Cellarium Cloud database to query by.
        :param metadata_feature_names: List of metadata feature names to include in the response.

        :return: List of cells with metadata.
        """
        request_data = {
            "cas_cell_ids": cell_ids,
            "metadata_feature_names": metadata_feature_names,
        }
        return self.post_json(endpoint=endpoints.QUERY_CELLS_BY_IDS, data=request_data)

    async def async_annotate_cell_type_summary_statistics_strategy(
        self, adata_bytes: t.ByteString, model_name: str, include_extended_output: bool
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Request Cellarium Cloud infrastructure to annotate an input anndata file using cell type count statistics

        Refer to API Docs: {api_url}/api/docs#/cell-analysis/annotate_api_cellarium_cas_annotate_post

        :param adata_bytes: Validated anndata file
        :param model_name: Name of the model to use.
        :param include_extended_output: Whether to include dev metadata in the response.

        :return: A list of dictionaries with annotations.
        """
        request_data = {
            "model_name": model_name,
            "include_extended_output": str(include_extended_output),
        }
        return await self.async_post(
            endpoints.ANNOTATE_CELL_TYPE_SUMMARY_STATS_STRATEGY, file=adata_bytes, data=request_data
        )

    async def async_annotate_cell_type_ontology_aware_strategy_anndata(
        self, adata_bytes: t.ByteString, model_name: str, prune_threshold: float, weighting_prefactor: float
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Request Cellarium Cloud infrastructure to annotate an input anndata file using ontology-aware strategy

        Refer to API Docs:
        {api_url}/api/docs#/cell-operations/annotate_ontology_aware_strategy_api_cellarium_cell_operations_annotate_ontology_aware_strategy_post

        :param adata_bytes: Validated anndata file
        :param model_name: Name of the model to use.
        :param prune_threshold: Whether to normalize consensus result.
        :param weighting_prefactor: Weighting prefactor for the ontology-aware strategy.

        :return: A list of dictionaries with annotations.
        """
        request_data = {
            "model_name": model_name,
            "prune_threshold": str(prune_threshold),
            "weighting_prefactor": str(weighting_prefactor),
        }
        return await self.async_post(
            endpoints.ANNOTATE_CELL_TYPE_ONTOLOGY_AWARE_STRATEGY, file=adata_bytes, data=request_data
        )

    async def async_nearest_neighbor_search(
        self, adata_bytes: t.ByteString, model_name: str
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Request Cellarium Cloud infrastructure to search for nearest neighbors in an input anndata file

        Refer tp API Docs:
        {api_url}/api/docs#/cell-analysis/nearest_neighbor_search_api_cellarium_cas_nearest_neighbor_search_post

        :param adata_bytes: Validated anndata file
        :param model_name: Name of the model to use.

        :return: A list of dictionaries with annotations (query_id and cas_cell_id, distance).
        """
        request_data = {
            "model_name": model_name,
        }
        return await self.async_post(endpoints.NEAREST_NEIGHBOR_SEARCH, file=adata_bytes, data=request_data)
