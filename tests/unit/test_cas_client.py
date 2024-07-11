import asyncio
import datetime
import typing as t

import aiohttp
import anndata
import numpy as np
import pandas as pd
import requests
import scipy.sparse as sp
from mockito import ANY, captor, mock, unstub, verify, when
from mockito.matchers import ArgumentCaptor

from cellarium.cas import constants
from cellarium.cas.client import CASClient
from tests.unit.test_utils import async_context_manager_decorator, async_return

TEST_TOKEN = "test_token"
TEST_URL = "https://cas-host.io"
TEST_SCHEMA = "schema1"

NP_RANDOM_STATE = np.random.RandomState(0)


class TestCasClient:
    def setup_method(self) -> None:
        # This tracks all async post mocks that are created where the key is the url and the value is the
        # session mock object. Used for verifying that the calls were made.
        self.async_post_mocks: t.Dict[str, aiohttp.ClientSession] = {}
        # Create a new event loop for each test since this is an async-heavy test and this avoid unexpected behavior
        loop = asyncio.new_event_loop()
        asyncio.set_event_loop(loop)

    def teardown_method(self) -> None:
        unstub()
        self.async_post_mocks = {}

    def test_initialize(self):
        self._mock_constructor_calls()

        cas_client = CASClient(api_token=TEST_TOKEN, api_url=TEST_URL)
        # Verify that the expected header values were sent.
        # E.g. make sure that the tokens all match and that the session and action ids all match.
        sent_tokens: t.Set[str] = set()
        sent_sessions: t.Set[str] = set()
        sent_actions: t.Set[str] = set()
        expected_calls = [
            f"{TEST_URL}/api/cellarium-general/validate-token",
            f"{TEST_URL}/api/cellarium-general/application-info",
            f"{TEST_URL}/api/cellarium-general/list-models",
            f"{TEST_URL}/api/cellarium-general/feature-schemas",
        ]

        for url in expected_calls:
            header_captor: ArgumentCaptor[t.Dict[str, any]] = captor()
            verify(requests).get(url=url, headers=header_captor)
            headers = header_captor.value
            if constants.Headers.authorization in headers:
                sent_tokens.add(headers[constants.Headers.authorization])
            if constants.Headers.client_session_id in headers:
                sent_sessions.add(headers[constants.Headers.client_session_id])
            if constants.Headers.client_action_id in headers:
                sent_actions.add(headers[constants.Headers.client_action_id])

        assert len(sent_tokens) == 1
        assert f"Bearer {TEST_TOKEN}" in sent_tokens

        assert len(sent_sessions) == 1
        assert str(cas_client.client_session_id) in sent_sessions

        assert len(sent_actions) == 1

    def test_annotate_matrix_cell_type_summary_statistics_strategy(self):
        num_cells = 10
        self._mock_constructor_calls()
        self._mock_annotate_matrix_cell_type_summary_statistics_strategy_calls(num_cells=num_cells)

        cas_client = CASClient(api_token=TEST_TOKEN, api_url=TEST_URL)

        response = cas_client.annotate_matrix_cell_type_summary_statistics_strategy(
            matrix=self._mock_anndata_matrix(num_cells=num_cells), chunk_size=100
        )

        assert len(response) == num_cells
        self._verify_headers(
            urls=[
                f"{TEST_URL}/api/cellarium-general/validate-token",
                f"{TEST_URL}/api/cellarium-general/quota",
                f"{TEST_URL}/api/cellarium-general/application-info",
                f"{TEST_URL}/api/cellarium-general/list-models",
                f"{TEST_URL}/api/cellarium-general/feature-schemas",
                f"{TEST_URL}/api/cellarium-general/feature-schema/{TEST_SCHEMA}",
            ],
            async_post_urls=[
                f"{TEST_URL}/api/cellarium-cell-operations/annotate-cell-type-summary-statistics-strategy"
            ],
            active_session_id=cas_client.client_session_id,
            num__expected_actions=2,  # one for initialization and one for the annotation
        )

    def test_annotate_matrix_cell_type_summary_statistics_strategy_with_chunking(self):
        num_cells = 100
        self._mock_constructor_calls()
        self._mock_annotate_matrix_cell_type_summary_statistics_strategy_calls(num_cells=num_cells)

        cas_client = CASClient(api_token=TEST_TOKEN, api_url=TEST_URL)

        # This should cause 10 chunks to be sent
        response = cas_client.annotate_matrix_cell_type_summary_statistics_strategy(
            matrix=self._mock_anndata_matrix(num_cells=num_cells), chunk_size=10
        )

        assert len(response) == num_cells
        self._verify_headers(
            urls=[
                f"{TEST_URL}/api/cellarium-general/validate-token",
                f"{TEST_URL}/api/cellarium-general/quota",
                f"{TEST_URL}/api/cellarium-general/application-info",
                f"{TEST_URL}/api/cellarium-general/list-models",
                f"{TEST_URL}/api/cellarium-general/feature-schemas",
                f"{TEST_URL}/api/cellarium-general/feature-schema/{TEST_SCHEMA}",
            ],
            async_post_urls=[
                f"{TEST_URL}/api/cellarium-cell-operations/annotate-cell-type-summary-statistics-strategy"
            ],
            active_session_id=cas_client.client_session_id,
            num__expected_actions=2,  # one for initialization and one for the annotation
        )

    def test_annotate_matrix_cell_type_summary_statistics_strategy_with_several_calls(self):
        num_cells = 100
        self._mock_constructor_calls()
        self._mock_annotate_matrix_cell_type_summary_statistics_strategy_calls(num_cells=num_cells)

        cas_client = CASClient(api_token=TEST_TOKEN, api_url=TEST_URL)

        # This should cause 10 chunks to be sent
        response1 = cas_client.annotate_matrix_cell_type_summary_statistics_strategy(
            matrix=self._mock_anndata_matrix(num_cells=num_cells), chunk_size=10
        )
        assert len(response1) == num_cells

        response2 = cas_client.annotate_matrix_cell_type_summary_statistics_strategy(
            matrix=self._mock_anndata_matrix(num_cells=num_cells), chunk_size=10
        )
        assert len(response2) == num_cells

        self._verify_headers(
            urls=[
                f"{TEST_URL}/api/cellarium-general/validate-token",
                f"{TEST_URL}/api/cellarium-general/quota",
                f"{TEST_URL}/api/cellarium-general/application-info",
                f"{TEST_URL}/api/cellarium-general/list-models",
                f"{TEST_URL}/api/cellarium-general/feature-schemas",
                f"{TEST_URL}/api/cellarium-general/feature-schema/{TEST_SCHEMA}",
            ],
            async_post_urls=[
                f"{TEST_URL}/api/cellarium-cell-operations/annotate-cell-type-summary-statistics-strategy"
            ],
            active_session_id=cas_client.client_session_id,
            num__expected_actions=3,  # one for initialization and one for *each* annotation call
        )

    def _mock_constructor_calls(self):
        """
        Mocks the calls made by the CASClient constructor
        """
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/validate-token",
            status_code=200,
            response_body={"username": "foo", "email": "foo@bar.com"},
        )
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/application-info",
            status_code=200,
            response_body={"application_version": "1.0.0", "default_feature_schema": "foo"},
        )
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/list-models",
            status_code=200,
            response_body=[
                {
                    "model_name": "human-pca-001",
                    "schema_name": TEST_SCHEMA,
                    "is_default_model": True,
                    "embedding_dimension": 512,
                    "description": "",
                }
            ],
        )
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/feature-schemas",
            status_code=200,
            response_body=[
                {"schema_name": TEST_SCHEMA},
            ],
        )
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/validate-client-version",
            status_code=200,
            response_body={"is_valid": True, "min_version": "1.4.0"},
            method="post",
        )

    def _mock_annotate_matrix_cell_type_summary_statistics_strategy_calls(
        self, num_cells: int = 3, num_features: int = 3
    ):
        """
        Mocks the calls made by the CASClient to do an annotation call with the summary statistics strategy
        """
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/feature-schema/{TEST_SCHEMA}",
            status_code=200,
            response_body=[f"field{i}" for i in range(num_features)],
        )
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/quota",
            status_code=200,
            response_body={
                "user_id": 0,
                "quota": 1000,
                "remaining_quota": 1000,
                "quota_reset_date": datetime.datetime.today() + 7 * datetime.timedelta(days=1),
            },
        )
        self._mock_async_post_response(
            url=f"{TEST_URL}/api/cellarium-cell-operations/annotate-cell-type-summary-statistics-strategy",
            status_code=200,
            response_body=[
                {
                    "query_cell_id": f"cell{i}",
                    "matches": [
                        {
                            "cell_type": "erythrocyte",
                            "cell_count": 100,
                            "min_distance": 12.0,
                            "p25_distance": 11.0,
                            "median_distance": 10.0,
                            "p75_distance": 9.0,
                            "max_distance": 13.0,
                        }
                    ],
                }
                for i in range(num_cells)
            ],
        )

    def _mock_anndata_matrix(self, num_features: int = 3, num_cells: int = 3) -> anndata.AnnData:
        d = NP_RANDOM_STATE.randint(0, 500, size=(num_cells, num_features))
        X = sp.csr_matrix(d)
        obs = pd.DataFrame(index=[f"cell{i}" for i in range(num_cells)])
        return anndata.AnnData(X=X, obs=obs, dtype=np.float32)

    def _mock_response(
        self,
        url: str,
        status_code: int,
        response_body: t.Union[dict, list],
        method: str = "get",
        post_data: t.Union[dict, list] = None,
    ):
        response = mock(aiohttp.ClientResponse)
        response.status_code = status_code
        when(response).json().thenReturn(response_body)

        if method == "get":
            when(requests).get(url=url, headers=ANY).thenReturn(response)
        elif method == "post":
            body = post_data if post_data is not None else ANY
            when(requests).post(url=url, headers=ANY, json=body).thenReturn(response)
        else:
            raise ValueError(f"Unsupported method: {method}")

    def _mock_async_post_response(
        self, url: str, status_code: int, response_body: t.Union[dict, list], post_data: t.Union[dict, list] = None
    ):
        # Mock response
        response = mock(aiohttp.ClientResponse)
        response.status = status_code
        when(response).json().thenReturn(async_return(response_body))
        response = async_context_manager_decorator(response)

        # Mock the session
        session = mock(aiohttp.ClientSession)
        body = post_data if post_data is not None else ANY
        session = async_context_manager_decorator(session)
        when(session).post(url, headers=ANY, data=body).thenReturn(response)

        # Mock the aiohttp.ClientSession constructor
        when(aiohttp).ClientSession(connector=ANY, timeout=ANY).thenReturn(session)
        self.async_post_mocks[url] = session

    def _verify_headers(
        self, urls: t.List[str], async_post_urls: t.List[str], active_session_id: str, num__expected_actions: int
    ):
        """
        Verify that the expected header values were sent.
        E.g. make sure that the tokens all match and that the session and action ids all match.
        """
        sent_tokens: t.Set[str] = set()
        sent_sessions: t.Set[str] = set()
        sent_actions: t.Set[str] = set()

        for url in urls:
            header_captor: ArgumentCaptor[t.Dict[str, any]] = captor()
            verify(requests, atleast=1).get(url=url, headers=header_captor)
            headers = header_captor.value
            if constants.Headers.authorization in headers:
                sent_tokens.add(headers[constants.Headers.authorization])
            if constants.Headers.client_session_id in headers:
                sent_sessions.add(headers[constants.Headers.client_session_id])
            if constants.Headers.client_action_id in headers:
                sent_actions.add(headers[constants.Headers.client_action_id])

        for url in async_post_urls or []:
            header_captor: ArgumentCaptor[t.Dict[str, any]] = captor()
            verify(self.async_post_mocks[url], atleast=1).post(url, headers=header_captor, data=ANY)
            headers = header_captor.value
            if constants.Headers.authorization in headers:
                sent_tokens.add(headers[constants.Headers.authorization])
            if constants.Headers.client_session_id in headers:
                sent_sessions.add(headers[constants.Headers.client_session_id])
            if constants.Headers.client_action_id in headers:
                sent_actions.add(headers[constants.Headers.client_action_id])

        assert len(sent_tokens) == 1
        assert f"Bearer {TEST_TOKEN}" in sent_tokens

        assert len(sent_sessions) == 1
        assert str(active_session_id) in sent_sessions

        assert len(sent_actions) == num__expected_actions
