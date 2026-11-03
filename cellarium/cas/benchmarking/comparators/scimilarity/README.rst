SCimilarity Benchmarking
=========================

This comparator runs official unconstrained ``model_v1.1`` with SCimilarity API
``0.4.1`` and produces the same four-file annotate directory used by
``cellarium-cas benchmark all``. The model and reference are mounted at runtime,
not baked into the image. The official workflow requires ``linux/amd64`` and at
least 64 GB RAM.

Directory structure::

    Dockerfile
    run.py                            # H5AD to ranked SCimilarity predictions
    helpers/
        map_scimilarity_to_cas_labels.py
        scimilarity_to_cas_annotation.py

Inputs and inference
--------------------

Input genes must be symbols. The runner uses ``feature_name``, then
``gene_symbol``, then ``var_names``; it takes raw counts from ``layers["counts"]``
or, only when that layer is absent, ``X``. Thus fallback ``X`` must be raw counts,
not normalized expression. Duplicate symbols are consolidated. Inputs with fewer
than 5,000 model genes fail.

The model remains unconstrained: it never calls ``safelist_celltypes`` or
``blocklist_celltypes``. Per-cell scores are normalized inverse-distance vote
shares across the 50 nearest neighbors, not model probabilities. The model group
is always ``scimilarity_model_v1.1``.

Build and run
-------------

Build on an amd64 Docker host::

    docker build --platform linux/amd64 \
        -t run-scimilarity:local \
        cellarium/cas/benchmarking/comparators/scimilarity

Download and unpack official ``model_v1.1`` under a local model directory, then
mount that directory and input data::

    docker run --rm --platform linux/amd64 \
        -v /path/to/models:/models:ro \
        -v /path/to/data:/data \
        run-scimilarity:local \
        --model-path /models/model_v1.1 \
        --input-h5ad /data/my_dataset.h5ad \
        --output-dir-template /data/scimilarity_{stem} \
        --top-k 5

Repeat ``--input-h5ad`` to retain the loaded reference across inputs. Each output
directory contains ``scimilarity_predictions.csv`` and ``scimilarity_stats.json``.
The prediction CSV has barcode index rows and interleaved
``predicted_label_k`` / ``predicted_score_k`` columns.

Crosswalk
---------

Create a SCimilarity-specific CSV with ``scimilarity_label``, ``cl_id``, and
``cl_label``. Do not reuse the scTab CSV wholesale: corpora and output vocabularies
differ. Exact validated source-label rows may be reused.

Extract source labels from the downloaded model::

    from scimilarity import CellAnnotation

    ca = CellAnnotation(model_path="/models/model_v1.1")
    for label in sorted(ca.classes):
        print(label)

Map each printed label once to an extant Cell Ontology primary term in the exact
benchmark ``ontology_resource.json``. Resolve unmatched labels through current
synonyms, obsolete-term ``replaced_by`` metadata, and human-specific variants;
record manual rationale outside this repository. Never fuzzy-match labels.

Before use, require every model class exactly once, unique source labels,
``cl_id`` values matching ``CL:\\d{7}``, every ID present in the ontology resource,
and each ``cl_label`` equal to that resource's primary label.

Convert and benchmark
---------------------

Run conversion from project root after creating the crosswalk::

    from cellarium.cas.benchmarking.comparators.scimilarity.helpers.scimilarity_to_cas_annotation import (
        scimilarity_to_cas_annotation,
    )

    scimilarity_to_cas_annotation(
        predictions_csv_path="/data/scimilarity_my_dataset/scimilarity_predictions.csv",
        h5ad_path="/data/my_dataset.h5ad",
        crosswalk_csv_path="/data/scimilarity_crosswalk.csv",
        crosswalk_label_col="scimilarity_label",
        crosswalk_cl_id_col="cl_id",
        crosswalk_cl_label_col="cl_label",
        ontology_resource_path="/data/ontology_resource.json",
        output_dir="/data/scimilarity_annotate_dir",
        top_k=5,
    )

The helper writes ``inferred_labels.csv``, ``metadata.json``,
``ontology_response.json``, and ``ontology_resource.json``. It rejects missing
barcodes, crosswalk labels, and ontology terms; duplicate Cell Ontology terms per
cell have scores summed and are reranked.

Place this directory beside CAS, scTab, or Azimuth annotate directories for the
same H5AD, then run::

    cellarium-cas benchmark all \
        --annotate-dirs /data \
        --output-dir /data/benchmark_results \
        --gt-label <ground_truth_obs_column> \
        --inferred-label cas_cell_type_name_1 \
        --f-measure-top-k 5
