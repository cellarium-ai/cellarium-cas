import asyncio
import datetime
import functools
import math
import operator
import pkgutil
import time
import typing as t
import warnings
from contextlib import contextmanager
from uuid import UUID, uuid4

import anndata
from anndata import ImplicitModificationWarning
from deprecated import deprecated

from cellarium.cas.logging import logger
from cellarium.cas.service import action_context_manager

from . import _io, constants, exceptions, models, preprocessing, service, settings, version


@contextmanager
def suppress_implicit_modification_warning():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ImplicitModificationWarning)
        yield


CHUNK_SIZE_ANNOTATE_DEFAULT = 1000
CHUNK_SIZE_SEARCH_DEFAULT = 500
DEFAULT_PRUNE_THRESHOLD = 0.05
DEFAULT_WEIGHTING_PREFACTOR = 1.0

FEEDBACK_TEMPLATE = pkgutil.get_data(__name__, "resources/feedback_template.html").decode("utf-8")


class CASClient:
    """
    Client that is designed to communicate with the Cellarium Cloud Backend.

    :param api_token: API token issued by the Cellarium team.
    :param api_url: URL of the Cellarium Cloud Backend. Should be left blank in most cases.
    :param num_attempts_per_chunk: Number of attempts the client should make to annotate each chunk. |br|
        `Default:` ``3``
    """

    def __print_models(self, models):
        s = "Allowed model list in Cellarium CAS:\n"
        for model in models:
            model_name = model["model_name"]
            description = model["description"]
            model_schema = model["schema_name"]
            embedding_dimension = model["embedding_dimension"]
            if model["is_default_model"]:
                model_name += " (default)"

            s += f"  - {model_name}\n    Description: {description}\n    Schema: {model_schema}\n    Embedding dimension: {embedding_dimension}\n"

        self.__print(s)

    @action_context_manager()
    def __init__(
        self,
        api_token: str,
        api_url: str = settings.CELLARIUM_CLOUD_BACKEND_URL,
        num_attempts_per_chunk: int = settings.NUM_ATTEMPTS_PER_CHUNK_DEFAULT,
    ) -> None:
        self.client_session_id: UUID = uuid4()
        self.cas_api_service = service.CASAPIService(
            api_token=api_token, api_url=api_url, client_session_id=self.client_session_id
        )

        logger.info(f"Connecting to the Cellarium Cloud backend with session {self.client_session_id}...")
        self.user_info = self.cas_api_service.validate_token()
        username = self.user_info["username"]
        self.__print(f"User is {username}")
        self.should_show_feedback = True
        if "should_ask_for_feedback" in self.user_info:
            self.should_show_feedback = self.user_info["should_ask_for_feedback"]

        self.validate_version()

        # Retrieving General Info
        application_info = self.cas_api_service.get_application_info()

        self.model_objects_list = self.cas_api_service.get_model_list()
        self.__allowed_models_list = [x["model_name"] for x in self.model_objects_list]
        self._model_name_obj_map = {x["model_name"]: x for x in self.model_objects_list}
        self.default_model_name = [x["model_name"] for x in self.model_objects_list if x["is_default_model"]][0]

        self.feature_schemas = self.cas_api_service.get_feature_schemas()
        self._feature_schemas_cache = {}

        self.num_attempts_per_chunk = num_attempts_per_chunk
        self.__print(f"Authenticated in Cellarium Cloud v. {application_info['application_version']}")

        self.__print_models(self.model_objects_list)

    @property
    def allowed_models_list(self):
        """
        List of models in Cellarium CAS that can be used to annotate.
        """
        return self.__allowed_models_list

    @staticmethod
    def __get_number_of_chunks(adata, chunk_size):
        return math.ceil(len(adata) / chunk_size)

    @staticmethod
    def __get_timestamp() -> str:
        return datetime.datetime.now().strftime("%H:%M:%S.%f")[:-3]

    def __print(self, str_to_print: str) -> None:
        print(f"* [{self.__get_timestamp()}] {str_to_print}")

    def __render_feedback_link(self):
        try:
            if settings.is_interactive_environment() and self.should_show_feedback:
                # only import IPython if we are in an interactive environment
                from IPython.display import HTML, display

                display(HTML(FEEDBACK_TEMPLATE.format(link=self.cas_api_service.get_feedback_answer_link())))
        except ModuleNotFoundError:
            pass

    def feedback_opt_out(self):
        self.should_show_feedback = False
        self.user_info = self.cas_api_service.feedback_opt_out()
        self.__print("Successfully opted out. You will no longer receive requests to provide feedback.")

    def validate_version(self):
        """
        Validate that this version of the client library is compatible with the selected server.
        """
        client_version = version.get_version()
        version_validation_info = self.cas_api_service.validate_version(version_str=client_version)
        if version_validation_info["is_valid"]:
            self.__print(f"Client version {client_version} is compatible with selected server.")
        else:
            raise exceptions.ClientTooOldError(
                f"Client version {client_version} is older than the minimum version for this server {version_validation_info['min_version']}. "
                f"Please update the client to the latest version using 'pip install cellarium-cas --upgrade'."
            )

    def validate_and_sanitize_input_data(
        self,
        adata: anndata.AnnData,
        cas_model_name: str,
        count_matrix_name: constants.CountMatrixInput,
        feature_ids_column_name: str,
        feature_names_column_name: t.Optional[str] = None,
    ) -> anndata.AnnData:
        """
        Validate and sanitize input :class:`anndata.AnnData` instance according to a specified feature schema associated
        with a particular model.

        :param adata: :class:`anndata.AnnData` instance to annotate
        :param cas_model_name: The model associated with the schema used for sanitizing. |br|
            `Allowed Values:` Model name from the :attr:`allowed_models_list` list keyword, which refers to the
            default selected model in the Cellarium backend. |br|
        :param count_matrix_name:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choice of either ``"X"``  or ``"raw.X"`` in order to use ``adata.X`` or ``adata.raw.X`` |br|
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
        :param feature_names_column_name: Column name where to obtain feature names (symbols).
            feature names wouldn't be mapped if value is ``None`` |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``None``
        :return: Validated and sanitized instance of :class:`anndata.AnnData`
        """
        cas_model_obj = self._model_name_obj_map[cas_model_name]
        feature_schema_name = cas_model_obj["schema_name"]
        try:
            cas_feature_schema_list = self._feature_schemas_cache[feature_schema_name]
        except KeyError:
            cas_feature_schema_list = self.cas_api_service.get_feature_schema_by(name=feature_schema_name)
            self._feature_schemas_cache[feature_schema_name] = cas_feature_schema_list

        # Make a copy of the input anndata object that we can sanitize
        new_adata = adata.copy()

        # Apply preprocessing steps that are necessary before validation and sanitization
        preprocessing.pre_sanitize(adata=new_adata, count_matrix_input=count_matrix_name)

        try:
            preprocessing.validate(
                adata=new_adata,
                cas_feature_schema_list=cas_feature_schema_list,
                feature_ids_column_name=feature_ids_column_name,
                count_matrix_input=count_matrix_name,
            )
        except exceptions.DataValidationError as e:
            if e.extra_features > 0:
                self.__print(
                    f"The input data matrix has {e.extra_features} extra features compared to '{feature_schema_name}' "
                    f"CAS schema ({len(cas_feature_schema_list)}). "
                    f"Extra input features will be dropped."
                )
            if e.missing_features > 0:
                self.__print(
                    f"The input data matrix has {e.missing_features} missing features compared to "
                    f"'{feature_schema_name}' CAS schema ({len(cas_feature_schema_list)}). "
                    f"Missing features will be imputed with zeros."
                )
            if e.extra_features == 0 and e.missing_features == 0:
                self.__print(
                    f"Input datafile has all the necessary features as {feature_schema_name}, but it's still "
                    f"incompatible because of the different order. The features will be reordered "
                    f"according to {feature_schema_name}..."
                )
                self.__print(
                    f"The input data matrix contains all of the features specified in '{feature_schema_name}' "
                    f"CAS schema but in a different order. The input features will be reordered according to "
                    f"'{feature_schema_name}'"
                )
            return preprocessing.sanitize(
                adata=new_adata,
                cas_feature_schema_list=cas_feature_schema_list,
                count_matrix_input=count_matrix_name,
                feature_ids_column_name=feature_ids_column_name,
                feature_names_column_name=feature_names_column_name,
            )
        else:
            self.__print(f"The input data matrix conforms with the '{feature_schema_name}' CAS schema.")
            return new_adata

    def print_user_quota(self) -> None:
        """
        Print the user's quota information
        """
        user_quota = self.cas_api_service.get_user_quota()
        lifetime_quota = user_quota["lifetime_quota"] or "unlimited"
        remaining_lifetime_quota = user_quota["remaining_lifetime_quota"] or "unlimited"
        self.__print(
            f"Weekly quota: {user_quota['weekly_quota']}, Remaining weekly quota: {user_quota['remaining_weekly_quota']}, "
            f"Weekly quota reset date: {user_quota['quota_reset_date']}\n"
            f"Lifetime quota: {lifetime_quota}, Remaining lifetime quota: {remaining_lifetime_quota} "
        )

    def __get_async_sharded_request_callback(
        self,
        results: t.List,
        chunk_index: int,
        chunk_start_i: int,
        chunk_end_i: int,
        semaphore: asyncio.Semaphore,
        service_request_callback: t.Callable,
    ) -> t.Callable:
        """
        A wrapper around POST request that handles HTTP and client errors, and resubmits the task if necessary.

        In case of HTTP 500, 503, 504 or ClientError print a message and resubmit the task.
        In case of HTTP 401 print a message to check the token.
        In case of any other error print a message.

        :param results: Results list that needs to be used to inplace the response from the server
        :param chunk_index: Consequent number of the chunk (e.g. Chunk 1, Chunk 2)
        :param chunk_start_i: Index pointing to the main adata file start position of the current chunk
        :param chunk_end_i: Index pointing to the main adata file end position of the current chunk
        :param semaphore: Semaphore object to limit the number of concurrent requests at a time
        :param service_request_callback: Callback function to execute (should be one of the async methods of the
            :class:`CASAPIService`)

        :return: A callback function that can be used to submit a request to the backend
        """

        async def sharded_request_task(**callback_kwargs):
            async with semaphore:
                retry_delay = settings.START_RETRY_DELAY

                for _ in range(self.num_attempts_per_chunk):
                    try:
                        results[chunk_index] = await service_request_callback(**callback_kwargs)

                    except (exceptions.HTTPError5XX, exceptions.HTTPClientError) as e:
                        self.__print(str(e))
                        self.__print(
                            f"Resubmitting chunk #{chunk_index + 1:2.0f} ({chunk_start_i:5.0f}, "
                            f"{chunk_end_i:5.0f}) to CAS ..."
                        )
                        await asyncio.sleep(retry_delay)
                        retry_delay = min(retry_delay * 2, settings.MAX_RETRY_DELAY)
                        continue
                    except exceptions.HTTPError401:
                        self.__print("Unauthorized token. Please check your API token or request a new one.")
                        break
                    except exceptions.HTTPError403 as e:
                        self.__print(str(e))
                        break
                    except Exception as e:
                        self.__print(f"Unexpected error: {e.__class__.__name__}; Message: {str(e)}")
                        break
                    else:
                        self.__print(
                            f"Received the result for cell chunk #{chunk_index + 1:2.0f} ({chunk_start_i:5.0f}, "
                            f"{chunk_end_i:5.0f}) ..."
                        )
                        break

        return sharded_request_task

    @suppress_implicit_modification_warning()
    def __async_sharded_request(
        self,
        adata: anndata.AnnData,
        chunk_size: int,
        request_callback: t.Callable,
        request_callback_kwargs: t.Dict[str, t.Any],
    ) -> t.List[t.Dict[str, t.Any]]:
        async def sharded_request():
            i, j = 0, chunk_size
            tasks = []
            semaphore = asyncio.Semaphore(settings.MAX_NUM_REQUESTS_AT_A_TIME)
            number_of_chunks = self.__get_number_of_chunks(adata, chunk_size=chunk_size)
            results = [[] for _ in range(number_of_chunks)]

            for chunk_index in range(number_of_chunks):
                chunk = adata[i:j, :]
                chunk_start_i = i
                chunk_end_i = i + len(chunk)
                self.__print(
                    f"Submitting cell chunk #{chunk_index + 1:2.0f} ({chunk_start_i:5.0f}, {chunk_end_i:5.0f}) "
                    f"to CAS ..."
                )

                chunk_bytes = _io.adata_to_bytes(adata=chunk)
                async_sharded_request_callback = self.__get_async_sharded_request_callback(
                    results=results,
                    chunk_index=chunk_index,
                    chunk_start_i=chunk_start_i,
                    chunk_end_i=chunk_end_i,
                    semaphore=semaphore,
                    service_request_callback=request_callback,
                )
                async_sharded_request_task = asyncio.create_task(
                    async_sharded_request_callback(adata_bytes=chunk_bytes, **request_callback_kwargs)
                )
                tasks.append(async_sharded_request_task)
                i = j
                j += chunk_size

            await asyncio.wait(tasks)
            return functools.reduce(operator.iconcat, results, [])

        return asyncio.run(sharded_request())

    def __postprocess_sharded_response(
        self, query_response: t.List[t.Dict[str, t.Any]], adata: anndata.AnnData, query_item_list_key: str
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Postprocess results by matching the order of cells in the response with the order of cells in the input

        :param query_response: List of dictionaries with annotations for each of the cells from input adata
        :param adata: :class:`anndata.AnnData` instance to annotate
        :param query_item_list_key: Key in the dictionary that contains the list of items (e.g. annotation matches,
        search result items)

        :return: A list of dictionaries with annotations for each of the cells from input adata, ordered by the input
        """

        processed_response = []
        query_response_hash = {x["query_cell_id"]: x for x in query_response}

        num_unannotated_cells = 0

        for query_cell_id in adata.obs.index:
            try:
                query_item = query_response_hash[query_cell_id]
            except KeyError:
                query_item = {"query_cell_id": query_cell_id, query_item_list_key: []}
                num_unannotated_cells += 1

            processed_response.append(query_item)

        if num_unannotated_cells > 0:
            self.__print(f"{num_unannotated_cells} cells were not processed by CAS")

        return processed_response

    def __postprocess_annotations(
        self, query_response: t.List[t.Dict[str, t.Any]], adata: anndata.AnnData
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Postprocess results by matching the order of cells in the response with the order of cells in the input

        :param query_response: List of dictionaries with annotations for each of the cells from input adata
        :param adata: :class:`anndata.AnnData` instance to annotate

        :return: A list of dictionaries with annotations for each of the cells from input adata
        """
        return self.__postprocess_sharded_response(
            query_response=query_response,
            adata=adata,
            query_item_list_key="matches",
        )

    def __postprocess_nearest_neighbor_search_response(
        self, query_response: t.List[t.Dict[str, t.Any]], adata: anndata.AnnData
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Postprocess nearest neighbor search response by matching the order of cells in the response with the order of
        cells in the input

        :param query_response: List of dictionaries with annotations for each of the cells from input adata
        :param adata: :class:`anndata.AnnData` instance to annotate

        :return: A list of dictionaries with nearest neighbor search results for each of the cells from input adata
        """
        return self.__postprocess_sharded_response(
            query_response=query_response,
            adata=adata,
            query_item_list_key="neighbors",
        )

    @staticmethod
    def __postprocess_query_cells_by_ids_response(
        query_response: t.List[t.Dict[str, t.Any]],
    ) -> t.List[t.Dict[str, t.Any]]:
        """
        Postprocess query cells by ids response by removing None values from cell metadata items in response

        :param query_response: List of dictionaries with annotations for each of the cells from input adata

        Returns:
        """
        processed_response = []

        for cell_metadata_item in query_response:
            cell_metadata_item_processed = {}

            for feature_name, value in cell_metadata_item.items():
                if value is None:
                    continue

                cell_metadata_item_processed[feature_name] = value

            processed_response.append(cell_metadata_item_processed)

        return processed_response

    def __validate_cells_under_quota(self, cell_count: int) -> None:
        """
        Validate the number of cells in the input data does not exceed remaining user quota

        :param cell_count: Number of cells in the input data
        """
        user_quota = self.cas_api_service.get_user_quota()
        if user_quota["remaining_lifetime_quota"] is not None and cell_count > user_quota["remaining_lifetime_quota"]:
            raise exceptions.QuotaExceededError(
                f"Number of cells in the input data ({cell_count}) exceeds the user's remaining lifetime quota ({user_quota['remaining_lifetime_quota']}).  "
                f"If you would like to discuss removing your lifetime quota, please reach out to cas-support@cellarium.ai"
            )
        elif cell_count > user_quota["remaining_weekly_quota"]:
            raise exceptions.QuotaExceededError(
                f"Number of cells in the input data ({cell_count}) exceeds the user's remaining quota ({user_quota['remaining_weekly_quota']}).  "
                f"The user's quota will be reset to {user_quota['weekly_quota']} on {user_quota['quota_reset_date']}."
            )

    def __prepare_input_for_sharded_request(
        self,
        adata: anndata.AnnData,
        count_matrix_input: constants.CountMatrixInput = constants.CountMatrixInput.X,
        feature_ids_column_name: str = "index",
        cas_model_name: t.Optional[str] = None,
        feature_names_column_name: t.Optional[str] = None,
    ) -> anndata.AnnData:
        """
        Prepare input data for sharded request. Validates model name and input data, and sanitizes the input data

        :param adata: :class:`anndata.AnnData` instance to annotate
        :param count_matrix_input:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values: Choices from enum :class:`constants.CountMatrixInput` |br|
            `Default:` ``"CountMatrixInput.X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :param cas_model_name: Model name to use for annotation. |br|
            `Allowed Values:` Model name from the :attr:`allowed_models_list` list or ``None``, which refers to the
            default selected model in the Cellarium backend. |br|
            `Default:` ``None``
        :param feature_names_column_name: Column name where to obtain feature names (symbols).
            feature names wouldn't be mapped if value is ``None`` |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``None``

        :return: A list of dictionaries with annotations for each of the cells from input adata
        """
        if cas_model_name not in self.allowed_models_list and cas_model_name is not None:
            raise ValueError(
                f"Model name '{cas_model_name}' is not in the list of allowed models. "
                f"Please use one of the following: {self.allowed_models_list} or `None` for the default model."
            )
        cas_model_name = self.default_model_name if cas_model_name is None else cas_model_name
        cas_model = self._model_name_obj_map[cas_model_name]
        cas_model_name = cas_model["model_name"]

        self.__print(f"Cellarium CAS (Model ID: {cas_model_name})")
        self.__print(f"Total number of input cells: {len(adata)}")

        self.__validate_cells_under_quota(cell_count=len(adata))

        return self.validate_and_sanitize_input_data(
            adata=adata,
            cas_model_name=cas_model_name,
            count_matrix_name=count_matrix_input,
            feature_ids_column_name=feature_ids_column_name,
            feature_names_column_name=feature_names_column_name,
        )

    @deprecated(version="1.4.3", reason="Use :meth:`annotate_matrix_cell_type_summary_statistics_strategy` instead")
    @action_context_manager()
    def annotate_anndata(
        self,
        adata: "anndata.AnnData",
        chunk_size=CHUNK_SIZE_ANNOTATE_DEFAULT,
        cas_model_name: str = "default",
        count_matrix_input: constants.CountMatrixInput = constants.CountMatrixInput.X,
        feature_ids_column_name: str = "index",
        feature_names_column_name: t.Optional[str] = None,
        include_dev_metadata: bool = False,
    ) -> models.CellTypeSummaryStatisticsResults:
        """
        Send an instance of :class:`anndata.AnnData` to the Cellarium Cloud backend for annotations. The function
        splits the ``adata`` into smaller chunks and asynchronously sends them to the backend API service. Each chunk is
        of equal size, except for the last one, which may be smaller. The backend processes these chunks in parallel.

        :param adata: :class:`anndata.AnnData` instance to annotate
        :param chunk_size: Size of chunks to split on
        :param cas_model_name: Model name to use for annotation. |br|
            `Allowed Values:` Model name from the :attr:`allowed_models_list` list or ``"default"``
            keyword, which refers to the default selected model in the Cellarium backend. |br|
            `Default:` ``"default"``
        :param count_matrix_input:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choices from enum :class:`~.CountMatrixInput` |br|
            `Default:` ``"CountMatrixInput.X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :param feature_names_column_name: Column name where to obtain feature names (symbols).
            feature names wouldn't be mapped if value is ``None`` |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``None``
        :param include_dev_metadata: Boolean indicating whether to include a breakdown of the number of cells
            by dataset

        :return: A :class:`~.models.CellTypeSummaryStatisticsResults` object with annotations for each of the cells from the
            adata input
        """
        cas_model_name = self.default_model_name if cas_model_name == "default" else cas_model_name

        start = time.time()
        adata = self.__prepare_input_for_sharded_request(
            adata=adata,
            cas_model_name=cas_model_name,
            count_matrix_input=count_matrix_input,
            feature_ids_column_name=feature_ids_column_name,
            feature_names_column_name=feature_names_column_name,
        )
        results = self.__async_sharded_request(
            adata=adata,
            chunk_size=chunk_size,
            request_callback=self.cas_api_service.async_annotate_cell_type_summary_statistics_strategy,
            request_callback_kwargs={
                "model_name": cas_model_name,
                "include_extended_output": include_dev_metadata,
            },
        )
        result = self.__postprocess_annotations(results, adata)
        # cast the object to the correct type
        result = models.CellTypeSummaryStatisticsResults(data=result)
        self.__print(f"Total wall clock time: {f'{time.time() - start:10.4f}'} seconds")
        self.__render_feedback_link()
        return result

    @deprecated(version="1.4.3", reason="Use :meth:`annotate_matrix_cell_type_statistics_strategy` instead")
    @action_context_manager()
    def annotate_anndata_file(
        self,
        filepath: str,
        chunk_size=CHUNK_SIZE_ANNOTATE_DEFAULT,
        cas_model_name: str = "default",
        count_matrix_input: constants.CountMatrixInput = constants.CountMatrixInput.X,
        feature_ids_column_name: str = "index",
        feature_names_column_name: t.Optional[str] = None,
        include_dev_metadata: bool = False,
    ) -> models.CellTypeSummaryStatisticsResults:
        """
        Read the 'h5ad' file into a :class:`anndata.AnnData` matrix and apply the :meth:`annotate_anndata` method to it.

        :param filepath: Filepath of the local :class:`anndata.AnnData` matrix
        :param chunk_size: Size of chunks to split on
        :param cas_model_name: Model name to use for annotation. |br|
            `Allowed Values:` Model name from the :attr:`allowed_models_list` list or ``"default"``
            keyword, which refers to the default selected model in the Cellarium backend. |br|
            `Default:` ``"default"``
        :param count_matrix_input:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choices from enum :class:`~.CountMatrixInput` |br|
            `Default:` ``"CountMatrixInput.X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :param feature_names_column_name: Column name where to obtain feature names (symbols).
            feature names wouldn't be mapped if value is ``None`` |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``None``
        :param include_dev_metadata: Boolean indicating whether to include a breakdown of the number of cells
            per dataset

        :return: A :class:`~.models.CellTypeSummaryStatisticsResults` object with annotations for each of the cells from
            the input adata
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            adata = anndata.read_h5ad(filename=filepath)

        return self.annotate_anndata(
            adata=adata,
            chunk_size=chunk_size,
            cas_model_name=cas_model_name,
            count_matrix_input=count_matrix_input,
            feature_ids_column_name=feature_ids_column_name,
            feature_names_column_name=feature_names_column_name,
            include_dev_metadata=include_dev_metadata,
        )

    @deprecated(version="1.4.3", reason="Use :meth:`annotate_matrix_with_cell_type_statistics_strategy` instead")
    @action_context_manager()
    def annotate_10x_h5_file(
        self,
        filepath: str,
        chunk_size: int = CHUNK_SIZE_ANNOTATE_DEFAULT,
        cas_model_name: str = "default",
        count_matrix_input: constants.CountMatrixInput = constants.CountMatrixInput.X,
        feature_ids_column_name: str = "index",
        feature_names_column_name: t.Optional[str] = None,
        include_dev_metadata: bool = False,
    ) -> models.CellTypeSummaryStatisticsResults:
        """
        Parse the 10x 'h5' matrix and apply the :meth:`annotate_anndata` method to it.

        :param filepath: Filepath of the local 'h5' matrix
        :param chunk_size: Size of chunks to split on
        :param cas_model_name: Model name to use for annotation. |br|
            `Allowed Values:` Model name from the :attr:`allowed_models_list` list or ``"default"``
            keyword, which refers to the default selected model in the Cellarium backend. |br|
            `Default:` ``"default"``
        :param count_matrix_input:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choices from enum :class:`~.CountMatrixInput` |br|
            `Default:` ``"CountMatrixInput.X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :param feature_names_column_name: Column name where to obtain feature names (symbols).
            feature names wouldn't be mapped if value is ``None`` |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``None``
        :param include_dev_metadata: Boolean indicating whether to include a breakdown of the number of cells by dataset

        :return: A :class:`~.models.CellTypeSummaryStatisticsResults` object with annotations for each of the cells from
            the input adata
        """
        adata = _io.read_10x_h5(filepath)

        return self.annotate_anndata(
            adata=adata,
            chunk_size=chunk_size,
            cas_model_name=cas_model_name,
            count_matrix_input=count_matrix_input,
            feature_ids_column_name=feature_ids_column_name,
            feature_names_column_name=feature_names_column_name,
            include_dev_metadata=include_dev_metadata,
        )

    @action_context_manager()
    def annotate_matrix_cell_type_summary_statistics_strategy(
        self,
        matrix: t.Union[str, anndata.AnnData],
        chunk_size=CHUNK_SIZE_ANNOTATE_DEFAULT,
        count_matrix_input: constants.CountMatrixInput = constants.CountMatrixInput.X,
        feature_ids_column_name: str = "index",
        include_extended_statistics: bool = True,
        cas_model_name: t.Optional[str] = None,
        feature_names_column_name: t.Optional[str] = None,
    ) -> models.CellTypeSummaryStatisticsResults:
        """
        Send an instance of :class:`anndata.AnnData` to the Cellarium Cloud backend for annotations. The function
        splits the ``adata`` into smaller chunks and asynchronously sends them to the backend API service. Each chunk is
        of equal size, except for the last one, which may be smaller. The backend processes these chunks in parallel.

        :param matrix: Either path to a file (must be either `.h5` or `.h5ad`) or an :class:`anndata.AnnData` instance
            to annotate
        :param chunk_size: Size of chunks to split on
        :param count_matrix_input:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choices from enum :class:`~.CountMatrixInput` |br|
            `Default:` ``"CountMatrixInput.X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :param include_extended_statistics: Boolean indicating whether to include a breakdown of the number of cells
            by dataset
        :param cas_model_name: Model name to use for annotation. |br|
            `Allowed Values:` Model name from the :attr:`allowed_models_list` list or ``None``, which refers to the
            default selected model in the Cellarium backend. |br|
            `Default:` ``None``
        :param feature_names_column_name: Column name where to obtain feature names (symbols).
            feature names wouldn't be mapped if value is ``None`` |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``None``

        :return: A :class:`~.models.CellTypeSummaryStatisticsResults` object with annotations for each of the cells from
            the input adata
        """
        if isinstance(matrix, str):
            matrix = _io.read_h5_or_h5ad(filename=matrix)

        cas_model_name = self.default_model_name if cas_model_name is None else cas_model_name

        start = time.time()
        matrix = self.__prepare_input_for_sharded_request(
            adata=matrix,
            cas_model_name=cas_model_name,
            count_matrix_input=count_matrix_input,
            feature_ids_column_name=feature_ids_column_name,
            feature_names_column_name=feature_names_column_name,
        )

        results = self.__async_sharded_request(
            adata=matrix,
            chunk_size=chunk_size,
            request_callback=self.cas_api_service.async_annotate_cell_type_summary_statistics_strategy,
            request_callback_kwargs={
                "model_name": cas_model_name,
                "include_extended_output": include_extended_statistics,
            },
        )
        result = self.__postprocess_annotations(query_response=results, adata=matrix)
        # cast the object to the correct type
        result = models.CellTypeSummaryStatisticsResults(data=result)
        self.__print(f"Total wall clock time: {f'{time.time() - start:10.4f}'} seconds")
        self.__render_feedback_link()
        return result

    @action_context_manager()
    def annotate_matrix_cell_type_ontology_aware_strategy(
        self,
        matrix: t.Union[str, anndata.AnnData],
        chunk_size=CHUNK_SIZE_ANNOTATE_DEFAULT,
        count_matrix_input: constants.CountMatrixInput = constants.CountMatrixInput.X,
        feature_ids_column_name: str = "index",
        cas_model_name: t.Optional[str] = None,
        feature_names_column_name: t.Optional[str] = None,
        prune_threshold: float = DEFAULT_PRUNE_THRESHOLD,
        weighting_prefactor: float = DEFAULT_WEIGHTING_PREFACTOR,
    ) -> models.CellTypeOntologyAwareResults:
        """
        Send an instance of :class:`anndata.AnnData` to the Cellarium Cloud backend for annotations using ontology
        aware strategy . The function splits the ``adata`` into smaller chunks and asynchronously sends them to the
        backend API service. Each chunk is of equal size, except for the last one, which may be smaller. The backend
        processes these chunks in parallel.

        :param matrix: Either path to a file (must be either `.h5` or `.h5ad`) or an :class:`anndata.AnnData` instance
            to annotate
        :param chunk_size: Size of chunks to split on
        :param count_matrix_input:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choices from enum :class:`~.CountMatrixInput` |br|
            `Default:` ``"CountMatrixInput.X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :param cas_model_name: Model name to use for annotation. |br|
            `Allowed Values:` Model name from the :attr:`allowed_models_list` list or ``None``, which refers to the
            default selected model in the Cellarium backend. |br|
            `Default:` ``None``
        :param feature_names_column_name: Column name where to obtain feature names (symbols).
            feature names wouldn't be mapped if value is ``None`` |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``None``
        :param prune_threshold: Threshold score for pruning the ontology graph in the output
        :param weighting_prefactor: Weighting prefactor for the weight calculation. A larger absolute value of the
            weighting_prefactor results in a steeper decay (weights drop off more quickly as distance increases),
            whereas a smaller absolute value results in a slower decay

        :return: A :class:`~.models.CellTypeOntologyAwareResults` object with annotations for each of the cells from
            the input adata
        """
        if isinstance(matrix, str):
            matrix = _io.read_h5_or_h5ad(filename=matrix)

        cas_model_name = self.default_model_name if cas_model_name is None else cas_model_name

        start = time.time()
        matrix = self.__prepare_input_for_sharded_request(
            adata=matrix,
            cas_model_name=cas_model_name,
            count_matrix_input=count_matrix_input,
            feature_ids_column_name=feature_ids_column_name,
            feature_names_column_name=feature_names_column_name,
        )
        results = self.__async_sharded_request(
            adata=matrix,
            chunk_size=chunk_size,
            request_callback=self.cas_api_service.async_annotate_cell_type_ontology_aware_strategy_anndata,
            request_callback_kwargs={
                "model_name": cas_model_name,
                "prune_threshold": prune_threshold,
                "weighting_prefactor": weighting_prefactor,
            },
        )
        result = self.__postprocess_annotations(results, matrix)
        # cast the object to the correct type
        result = models.CellTypeOntologyAwareResults(data=result)
        self.__print(f"Total wall clock time: {f'{time.time() - start:10.4f}'} seconds")
        self.__render_feedback_link()
        return result

    @deprecated(version="1.4.3", reason="Use :meth:`search_matrix` instead")
    @action_context_manager()
    def search_anndata(
        self,
        adata: anndata.AnnData,
        chunk_size=CHUNK_SIZE_SEARCH_DEFAULT,
        cas_model_name: str = "default",
        count_matrix_input: constants.CountMatrixInput = constants.CountMatrixInput.X,
        feature_ids_column_name: str = "index",
        feature_names_column_name: t.Optional[str] = None,
    ) -> models.MatrixQueryResults:
        """
        Send an instance of :class:`anndata.AnnData` to the Cellarium Cloud backend for nearest neighbor search. The
        function splits the ``adata`` into smaller chunks and asynchronously sends them to the backend API service.
        Each chunk is of equal size, except for the last one, which may be smaller. The backend processes
        these chunks in parallel.

        :param adata: :class:`anndata.AnnData` instance to annotate
        :param chunk_size: Size of chunks to split on
        :param cas_model_name: Model name to use for annotation. |br|
            `Allowed Values:` Model name from the :attr:`allowed_models_list` list or ``"default"``
            keyword, which refers to the default selected model in the Cellarium backend. |br|
            `Default:` ``"default"``
        :param count_matrix_input:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choices from enum :class:`~.CountMatrixInput` |br|
            `Default:` ``"CountMatrixInput.X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :param feature_names_column_name: Column name where to obtain feature names (symbols).
            feature names wouldn't be mapped if value is ``None`` |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``None``

        :return: A :class:`~.models.MatrixQueryResults` object with search results for each of the cells from
            the input adata
        """
        if chunk_size > settings.MAX_CHUNK_SIZE_SEARCH_METHOD:
            raise ValueError("Chunk size greater than 500 not supported yet.")

        cas_model_name = self.default_model_name if cas_model_name == "default" else cas_model_name

        start = time.time()
        adata = self.__prepare_input_for_sharded_request(
            adata=adata,
            cas_model_name=cas_model_name,
            count_matrix_input=count_matrix_input,
            feature_ids_column_name=feature_ids_column_name,
            feature_names_column_name=feature_names_column_name,
        )
        results = self.__async_sharded_request(
            adata=adata,
            chunk_size=chunk_size,
            request_callback=self.cas_api_service.async_nearest_neighbor_search,
            request_callback_kwargs={"model_name": cas_model_name},
        )
        result = self.__postprocess_nearest_neighbor_search_response(results, adata)
        # cast the object to the correct type
        result = models.MatrixQueryResults(data=result)
        self.__print(f"Total wall clock time: {f'{time.time() - start:10.4f}'} seconds")
        return result

    @deprecated(version="1.4.3", reason="Use :meth:`search_matrix` instead")
    @action_context_manager()
    def search_10x_h5_file(
        self,
        filepath: str,
        chunk_size: int = CHUNK_SIZE_SEARCH_DEFAULT,
        cas_model_name: str = "default",
        count_matrix_input: constants.CountMatrixInput = constants.CountMatrixInput.X,
        feature_ids_column_name: str = "index",
        feature_names_column_name: t.Optional[str] = None,
    ) -> models.MatrixQueryResults:
        """
        Parse the 10x 'h5' matrix and apply the :meth:`search_anndata` method to it.

        :param filepath: Filepath of the local 'h5' matrix
        :param chunk_size: Size of chunks to split on
        :param cas_model_name: Model name to use for annotation. |br|
            `Allowed Values:` Model name from the :attr:`allowed_models_list` list or ``"default"``
            keyword, which refers to the default selected model in the Cellarium backend. |br|
            `Default:` ``"default"``
        :param count_matrix_input:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choices from enum :class:`~.CountMatrixInput` |br|
            `Default:` ``"CountMatrixInput.X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :param feature_names_column_name: Column name where to obtain feature names (symbols).
            feature names wouldn't be mapped if value is ``None`` |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``None``

        :return: A :class:`~.models.MatrixQueryResults` object with search results for each of the cells from
            the input adata
        """
        adata = _io.read_10x_h5(filepath)

        return self.search_anndata(
            adata=adata,
            chunk_size=chunk_size,
            cas_model_name=cas_model_name,
            count_matrix_input=count_matrix_input,
            feature_ids_column_name=feature_ids_column_name,
            feature_names_column_name=feature_names_column_name,
        )

    @action_context_manager()
    def search_matrix(
        self,
        matrix: t.Union[str, anndata.AnnData],
        chunk_size: int = CHUNK_SIZE_SEARCH_DEFAULT,
        count_matrix_input: constants.CountMatrixInput = constants.CountMatrixInput.X,
        feature_ids_column_name: str = "index",
        cas_model_name: t.Optional[str] = None,
        feature_names_column_name: t.Optional[str] = None,
    ) -> models.MatrixQueryResults:
        """
        Send an instance of :class:`anndata.AnnData` to the Cellarium Cloud backend for nearest neighbor search. The
        function splits the ``adata`` into smaller chunks and asynchronously sends them to the backend API service.
        Each chunk is of equal size, except for the last one, which may be smaller. The backend processes
        these chunks in parallel.

        :param matrix: Either path to a file (must be either `.h5` or `.h5ad`) or an :class:`anndata.AnnData` instance
            to annotate
        :param chunk_size: Size of chunks to split on
        :param count_matrix_input:  Where to obtain a feature expression count matrix from. |br|
            `Allowed Values:` Choices from enum :class:`~.CountMatrixInput` |br|
            `Default:` ``"CountMatrixInput.X"``
        :param feature_ids_column_name: Column name where to obtain Ensembl feature ids. |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``"index"``
        :param cas_model_name: Model name to use for annotation. |br|
            `Allowed Values:` Model name from the :attr:`allowed_models_list` list or ``None``, which refers to the
            default selected model in the Cellarium backend. |br|
            `Default:` ``None``
        :param feature_names_column_name: Column name where to obtain feature names (symbols).
            feature names wouldn't be mapped if value is ``None`` |br|
            `Allowed Values:` A value from ``adata.var.columns`` or ``"index"`` keyword, which refers to index
            column. |br|
            `Default:` ``None``

        :return: A :class:`~.models.MatrixQueryResults` object with search results for each of the cells from
            the input adata
        """
        if isinstance(matrix, str):
            matrix = _io.read_h5_or_h5ad(filename=matrix)

        if chunk_size > settings.MAX_CHUNK_SIZE_SEARCH_METHOD:
            raise ValueError(f"Chunk size greater than `{settings.MAX_CHUNK_SIZE_SEARCH_METHOD}` not supported yet.")

        cas_model_name = self.default_model_name if cas_model_name is None else cas_model_name

        start = time.time()
        matrix = self.__prepare_input_for_sharded_request(
            adata=matrix,
            cas_model_name=cas_model_name,
            count_matrix_input=count_matrix_input,
            feature_ids_column_name=feature_ids_column_name,
            feature_names_column_name=feature_names_column_name,
        )
        results = self.__async_sharded_request(
            adata=matrix,
            chunk_size=chunk_size,
            request_callback=self.cas_api_service.async_nearest_neighbor_search,
            request_callback_kwargs={"model_name": cas_model_name},
        )
        result = self.__postprocess_nearest_neighbor_search_response(results, matrix)
        # cast the object to the correct type
        result = models.MatrixQueryResults(data=result)
        self.__print(f"Total wall clock time: {f'{time.time() - start:10.4f}'} seconds")
        self.__render_feedback_link()
        return result

    @action_context_manager()
    def query_cells_by_ids(
        self, cell_ids: t.List[int], metadata_feature_names: t.List[constants.CellMetadataFeatures] = None
    ) -> models.CellQueryResults:
        """
        Query cells by their ids from a single anndata file with Cellarium CAS. Input file should be validated and
        sanitized according to the model schema.

        :param cell_ids: List of cell ids to query
        :param metadata_feature_names: List of metadata features to include in the response. |br|

        :return: A :class:`~.models.CellQueryResults` object with cell query results
        """
        results = self.cas_api_service.query_cells_by_ids(
            cell_ids=cell_ids,
            metadata_feature_names=metadata_feature_names,
        )
        result = self.__postprocess_query_cells_by_ids_response(query_response=results)
        # cast the object to the correct type
        result = models.CellQueryResults(data=result)
        self.__render_feedback_link()
        return result

    def validate_model_name(self, model_name: t.Optional[str] = None) -> None:
        """
        Validate if the model name provided is valid

        :param model_name: Model name to check
        :raises: ValueError if model name is not valid
        """
        if model_name not in self.allowed_models_list and model_name is not None:
            raise ValueError(
                f"Model name '{model_name}' is not in the list of allowed models. "
                f"Please use one of the following: {self.allowed_models_list} or `None` for the default model."
            )
