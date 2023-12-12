Changelog
#########

All notable changes to Cellarium CAS client will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.


1.4.1.dev - 2023-12-08
----------------------

Added
~~~~~
- Include kNN search method (#49)
- Include get cells by IDs method (#49)
- Include helper methods for visualization and demo
- Add model name validation method to :class:`clients.CASClient`
- Add sync POST method (using requests) to :class:`services.CASAPIService`
- Add `CHANGELOG.rst` file

Changed
~~~~~~~
- Reorganize :class:`CASClient` methods: factor out sharding logic
- Update `MAX_NUM_REQUESTS_AT_A_TIME` to 10
- Update default `chunk_size` in :meth:`annotate` methods to 1000
- Make :meth:`__validate_and_sanitize_input_data` method public (now it's a :meth:`validate_and_sanitize_input_data`) in CASClient
- Update backend API url to `https://cas-api-1-4-1-dev-vi7nxpvk7a-uc.a.run.app`
- Update `pyproject.toml` file to include scanpy optional dependencies

Removed
~~~~~~~
- Remove docs generation from CI/CD pipeline

File Structure Changes
~~~~~~~~~~~~~~~~~~~~~~
- Add `CHANGELOG.rst` file
- Add `requirements/scanpy.txt` file (optional requirements for scanpy related demos)
- Add `cellarium/cas/scanpy_utils.py` (Not necessary for the client methods, but useful for the demo)
- Remove `.github/actions/docs` folder (docs are now hosted on readthedocs)

Notes
~~~~~
- Users will need a new API token to use this version


