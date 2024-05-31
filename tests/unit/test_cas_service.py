import aiohttp
import pytest
import requests
from mockito import mock, unstub, when

from cellarium.cas import exceptions
from cellarium.cas.service import CASAPIService

TEST_TOKEN = "test_token"
TEST_URL = "https://cas-host.io"


class TestCasService:
    def setup_method(self) -> None:
        self.cas_service = CASAPIService(api_token=TEST_TOKEN, api_url=TEST_URL)

    def teardown_method(self) -> None:
        unstub()

    def test_validate_token(self):
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/validate-token",
            token=TEST_TOKEN,
            status_code=200,
            response_body={"username": "foo", "email": "foo@bar.com"},
        )

        self.cas_service.validate_token()

    def test_validate_token_invalid(self):
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/validate-token",
            token=TEST_TOKEN,
            status_code=401,
            response_body={"detail": "Invalid Token"},
        )

        with pytest.raises(exceptions.HTTPError401, match="Invalid Token"):
            self.cas_service.validate_token()

    def test_get_application_info(self):
        response_body = {"application_version": "1.0.0", "default_feature_schema": "foo"}
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/application-info",
            token=TEST_TOKEN,
            status_code=200,
            response_body=response_body,
        )

        assert self.cas_service.get_application_info() == response_body

    def test_get_feature_schemas(self):
        response_body = [
            {"schema_name": "schema1"},
            {"schema_name": "schema2"},
        ]
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/feature-schemas",
            token=TEST_TOKEN,
            status_code=200,
            response_body=response_body,
        )

        assert self.cas_service.get_feature_schemas() == ["schema1", "schema2"]

    def test_get_feature_schema_by(self):
        response_body = ["field1", "field2"]
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/feature-schema/schema1",
            token=TEST_TOKEN,
            status_code=200,
            response_body=response_body,
        )

        assert self.cas_service.get_feature_schema_by(name="schema1") == ["field1", "field2"]

    def test_get_feature_schema_by_not_found(self):
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/feature-schema/schema1",
            token=TEST_TOKEN,
            status_code=404,
            response_body={"detail": "Not found"},
        )

        with pytest.raises(exceptions.HTTPError404, match="Not found"):
            self.cas_service.get_feature_schema_by(name="schema1")

    def test_get_feature_schema_by_unauthorized(self):
        self._mock_response(
            url=f"{TEST_URL}/api/cellarium-general/feature-schema/schema1",
            token=TEST_TOKEN,
            status_code=403,
            response_body={"detail": "Unauthorized"},
        )

        with pytest.raises(exceptions.HTTPError403, match="Unauthorized"):
            self.cas_service.get_feature_schema_by(name="schema1")

    def _mock_response(self, url: str, token: str, status_code: int, response_body: dict | list):
        response = mock(aiohttp.ClientResponse)
        response.status_code = status_code
        when(response).json().thenReturn(response_body)
        when(requests).get(url=url, headers={"Authorization": f"Bearer {token}"}).thenReturn(response)
