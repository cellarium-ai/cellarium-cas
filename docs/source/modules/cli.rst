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
to also write ``inferred_labels.csv`` with per-cell top-k cell type assignments.
Pass ``--cluster-label`` to additionally compute cluster-level label calls.

.. code-block:: bash

    cellarium-cas annotate \
        --input-path cells.h5ad \
        --output-dir ./cas_output \
        --cluster-label leiden

.. click:: cellarium.cas.cli.annotate:annotate_command
   :prog: cellarium-cas annotate
   :nested: full

----

benchmark
---------

Evaluate annotation quality against a labelled reference dataset. Both subcommands
consume output directories produced by ``cellarium-cas annotate``.

.. click:: cellarium.cas.cli.benchmark:benchmark_group
   :prog: cellarium-cas benchmark
   :nested: full
