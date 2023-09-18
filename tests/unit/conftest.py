"""
This module contains configuration for pytest, including fixtures, hooks, plugins,
and other configurations that are shared across multiple test files.

`pytest` automatically discovers this file and uses it as a configuration source.

Attributes:
    Any global variables, fixtures, or configurations you define.

Functions:
    pytest_addoption(parser): Adds custom command line options.
    pytest_configure(config): Sets custom configurations using the provided options.

Fixtures:
    test_api_token(request): Provides access to the API token for testing.

Notes:
    Ensure that this file resides in a location where `pytest` can discover it,
    usually in the root of the test directory or subdirectories.

"""

import pytest


def pytest_addoption(parser: pytest.Parser):
    """Add `test_api_token` option to pytest args"""
    parser.addoption("--test_api_token", action="store", default="default_token", help="Token for API calls")


@pytest.fixture
def test_api_token(request: pytest.FixtureRequest) -> str:
    """Get `test_api_token` from parser's options"""
    return request.config.getoption("--test_api_token")
