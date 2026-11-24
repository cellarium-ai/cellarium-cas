:tocdepth: 3

Command-Line Interface
######################

``cellarium-cas`` ships a command-line interface for running annotation and benchmarking jobs
without writing Python. The API token is read from the ``CAS_API_TOKEN`` environment variable
by default, or passed explicitly via ``--cas-api-token``.

.. code-block:: bash

    export CAS_API_TOKEN="your_token_here"

----

annotate
--------

Annotate a local ``.h5ad`` file using the CAS ontology-aware strategy and write outputs
to a directory. The output directory will contain ``ontology_response.json``,
``ontology_resource.json``, and ``metadata.json`` by default. Pass ``--infer-labels``
to also write ``inferred_labels.csv`` with per-cell cell type assignments at four granularity positions.
Pass ``--cluster-label`` to additionally compute cluster-level label calls. Pass
``--knn-smoothing K`` to smooth CAS scores for each cell over its ``K`` nearest query-cell
neighbors while preserving single-cell resolution; this writes ``cas_annotation_knn_*`` columns.
Call selection uses ``compute_annotation_levels_*`` and the curated L1-L4 granularity lookup. The lookup is fetched from the
backend with the same ontology resource name as the cell ontology resource and cached on the
``CASClient``.

.. code-block:: bash

    cellarium-cas annotate \
        --input-path cells.h5ad \
        --output-dir ./cas_output \
        --cluster-label leiden \
        --knn-smoothing 30

.. click:: cellarium.cas.cli.annotate:annotate_command
   :prog: cellarium-cas annotate
   :nested: full

----

benchmark
---------

Evaluate annotation quality against a labelled reference dataset.  The benchmark
subcommands consume output directories produced by ``cellarium-cas annotate``.
All benchmark artifacts (confusion matrices, metric CSVs) are written to a single
``--output-dir`` that acts as the benchmarking workspace.  Annotate directories
are never modified.

The pipeline has four ordered steps plus a convenience command that runs all of them:

.. code-block:: bash

    # All-in-one
    cellarium-cas benchmark all \
        --annotate-dirs ./annotate_outputs \
        --output-dir    ./benchmark_results \
        --gt-label      cell_type_ontology_term_id \
        --inferred-label cas_annotation_name_1

    # Or step by step
    cellarium-cas benchmark confusion-matrix ...
    cellarium-cas benchmark aggregate        ...
    cellarium-cas benchmark f-measure        ...
    cellarium-cas benchmark hierarchical     ...

Keep ``--f-measure-top-k 1`` for ``cas_annotation_*`` columns: their numeric suffixes
are granularity levels, not alternative prediction ranks.

See :doc:`benchmarking` for full documentation of the pipeline and output columns.

.. click:: cellarium.cas.cli.benchmark:benchmark_group
   :prog: cellarium-cas benchmark
   :nested: full
