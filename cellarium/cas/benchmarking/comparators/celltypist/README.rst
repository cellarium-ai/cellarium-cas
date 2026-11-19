CellTypist Benchmarking
=======================

Raw per-cell CellTypist (``mode='best match'``, no majority voting) with
``celltypist==1.7.1``. A trusted local model ``.pkl`` is required per run; none is
bundled or downloaded (loading a model unpickles it).

Directory structure::

    Dockerfile
    run.py                            # H5AD to ranked native CellTypist predictions
    helpers/
        celltypist_to_cas_annotation.py

Inference
---------

Inputs need gene symbols (``feature_name``, ``gene_symbol``, then ``var_names``) and
raw counts (``layers["counts"]``, else ``X``). The runner sums duplicate symbols,
scales to 10,000 counts per cell, applies ``log1p``, and runs CellTypist. Ranks
are the native classes ordered by decision score; scores are the original sigmoid
probabilities (independent per class, not a distribution, and may saturate at 1)::

    docker build -t run-celltypist:local cellarium/cas/benchmarking/comparators/celltypist

    docker run --rm \
        -v /path/to/models:/models:ro -v /path/to/data:/data \
        run-celltypist:local \
        --model-path /models/Immune_All_Low.pkl \
        --input-h5ad /data/my_dataset.h5ad \
        --output-dir-template /data/celltypist_{stem} --top-k 5

Each output directory holds ``predicted_labels.csv``, ``celltypist_predictions.csv``
and ``celltypist_stats.json``. The benchmark group name is
``celltypist_<model-stem>_<model-sha256>_raw``.

Conversion
----------

``celltypist_to_cas_annotation`` builds the four-file CAS annotation directory. The
crosswalk must be model-specific with exactly one CL ID per native label; every
predicted label must map and every CL ID must exist in the ontology resource.
Labels sharing a CL ID within the exported top-k merge by maximum score. Run
``cellarium-cas benchmark all --inferred-label cas_cell_type_name_1`` on the
converted directories.
