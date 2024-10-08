Changelog
#########

All notable changes to Cellarium CAS client will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.


..
  The text in this block is where pre-release changes should live.
  On release, a commit should be created to copy the block just below with the new version number and date. 
  Then a new block should be created here for the next version.

  <pre-version> - <pre-date>
  --------------------------

  Changed
  ~~~~~~~

  Fixed
  ~~~~~

1.6.0 - 2024-10-08
------------------

Changed
~~~~~~~
- Added ability to enforce lifetime quotas
- Dependency vulnerability checks are performed by unit testing github action
- Improved handling of different root nodes in the circular tree plot visualization Dash App

Fixed
~~~~~
- Logging now uses Python logging 

1.5.1 - 2024-09-09
------------------

Changed
~~~~~~~
- Added ability to render documentation locally
- PR tests will now fail when doc is invalid
- Fixed bugs in documentation navigation
- Update Cell Type Ontology to the version used in cellxgene schema v5
- Explicitly typed all objects returned from the client

Fixed
~~~~~
- Added a pre-sanitization step to fix an issue where non-float32 anndata matrices would be submitted to CAS and fail
- Removed unused method parameters


1.4.13 - 2024-09-06
-------------------
**Warning: This version was compromised due to a wrong commit being tagged. It is not a proper release version. Please use version 1.5.1 or later.**


1.4.12 - 2024-08-01
-------------------

Changed
~~~~~~~
- Removed references to the development environment in the client code
- Updated client documentation
- Updated the demo notebook's documentation

1.4.11 - 2024-07-24
-------------------

Added
~~~~~
- Added docstrings in visualization and exposed it to sphinx docs


1.4.10 - 2024-07-23
-------------------

Fixed
~~~~~
- Bug fix in the circular tree plot visualization Dash App


1.4.9 - 2024-07-22
------------------

Added
~~~~~
- Added circular tree plot visualization Dash App

Changed
~~~~~~~
- Renamed ``data_preparation`` to ``preprocessing```
- Moved all preprocessing-related code to ``preprocessing`` submodule
- Added a ``postprocessing`` submodule
- Added a ``visualization`` submodule
- Renamed ``scanpy`` optional dependencies to ``vis`` (for all visualization-related dependencies)


1.4.8 - 2024-07-22
------------------

Changed
~~~~~~~

- Decrease `MAX_NUM_REQUESTS_AT_A_TIME` to 8


1.4.7 - 2024-07-09
------------------

Added
~~~~~
- Add check in client initialization to ensure the current version of the client code is compatbile with the selected CAS server

1.4.6 - 2024-06-12
------------------

Added
~~~~~
- Add text requesting feedback at the end of calls with instructions on opting out of feedback requests


1.4.5 - 2024-06-03
------------------

Added
~~~~~
- Add optional parameter to the client constructor to set the CAS API URL

1.4.4 - 2024-05-02
------------------

Added
~~~~~
- Add :meth:`~.CASClient.annotate_matrix_cell_type_summary_statistics_strategy` method to :class:`~.CASClient`
- Add :meth:`~.CASClient.annotate_matrix_cell_type_ontology_aware_strategy` method to :class:`~.CASClient`

Changed
~~~~~~~
- Deprecate :meth:`~.CASClient.annotate_anndata`, :meth:`~.CASClient.annotate_anndata_file`, :meth:`~.CASClient.annotate_10x_h5_file`, :meth:`~.CASClient.search_anndata`, and :meth:`~.CASClient.search_10x_h5_file`,  methods in :class:`~.CASClient`

File Structure Changes
~~~~~~~~~~~~~~~~~~~~~~
- No changes in file structure

1.4.3 - 2024-03-18
------------------

Added
~~~~~
- Fix total mrna umis for normalized data

Changed
~~~~~~~
- Handle different matrix types in the data preparation callbacks
- Update unit tests for the data preparation callbacks

1.4.2 - 2024-03-12
------------------

Changed
~~~~~~~
- Increase client HTTP request timeouts

1.4.1 - 2024-02-15
------------------

Added
~~~~~
- Include kNN search method (#49)
- Include get cells by IDs method (#49)
- Include helper methods for visualization and demo
- Add model name validation method to :class:`~.CASClient`
- Add sync POST method (using requests) to CASAPIService
- Add ``CHANGELOG.rst`` file
- Add settings module that chooses the correct settings file based on the environment according to current git version. Since now package will use development settings if it's tagged as a pre-release (alpha, beta, or release candidate (rc)), and production settings otherwise.
- Add version determination based on git tags
- Add callback methods to data_preparation module. Include total total_mrna_umis calculation as a callback before data sanitization

Changed
~~~~~~~
- Reorganize :class:`~.CASClient` methods: factor out sharding logic
- Update ``MAX_NUM_REQUESTS_AT_A_TIME`` to 25
- Update default ``chunk_size`` in :meth:`~.CASClient.annotate_anndata` methods to 1000
- Make ``__validate_and_sanitize_input_data`` method public (now it is :meth:`~.CASClient.validate_and_sanitize_input_data`) in :class:`~.CASClient`
- Update backend API url to point to the new API endpoints depending on the environment
- Update ``pyproject.toml`` file to include scanpy optional dependencies
- Restructure data_preparation into a module

Removed
~~~~~~~
- Remove docs generation from CI/CD pipeline

File Structure Changes
~~~~~~~~~~~~~~~~~~~~~~
- Add ``CHANGELOG.rst`` file
- Add ``requirements/scanpy.txt`` file (optional requirements for scanpy related demos)
- Add ``cellarium/cas/scanpy_utils.py`` (Not necessary for the client methods, but useful for the demo)
- Add ``cellarium/cas/settings`` directory, including ``__init__.py``, ``base.py``, ``development.py``, and ``production.py`` files
- Add cas/version.py file
- Add ``cellarium/cas/data_preparation`` directory, including ``__init__.py``, ``callbacks.py``, ``sanitizer.py`` and ``validator.py``` files
- Add ``tests/unit/test_data_preparation_callbacks.py`` file
- Add ``cellarium/cas/constants.py`` file
- Remove ``.github/actions/docs`` folder (docs are now hosted on readthedocs)

Notes
~~~~~
- Users will need a new API token to use this version
