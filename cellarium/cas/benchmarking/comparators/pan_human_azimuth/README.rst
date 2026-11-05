Pan-Human Azimuth Annotation
============================

This comparator runs official Pan-Human Azimuth model ``v1`` through
``panhumanpy==1.0.0`` (``AzimuthNN`` with refinement), writes its native fine,
medium, and broad annotations, and converts them to the CAS annotate-directory
format given a user-supplied crosswalk. It does not create the crosswalk. The
image is CPU-only TensorFlow and is built and supported on ``linux/amd64`` only.
It runs as nonroot ``appuser`` (UID 1000), so mounted output and cache
directories must be writable by UID 1000. No RAM requirement is established.

Model weights are not baked into the image. The first run downloads them over the
network to ``/home/appuser/.cache/panhumanpy``; mount a persistent host directory
there to avoid repeat downloads.

Directory structure::

    Dockerfile
    run.py                            # H5AD to native Pan-Human Azimuth annotations
    helpers/
        map_pan_human_azimuth_to_cas_labels.py
        pan_human_azimuth_to_cas_annotation.py

Inputs and inference
--------------------

Input genes must be symbols. The runner takes counts from ``layers["counts"]`` or,
only when that layer is absent, ``X``; ``.raw`` is never used. Counts must be finite,
nonnegative, and integer-valued (float storage of integers is fine); normalized or
fractional data fails. Symbols come from ``var["feature_name"]``, then
``var["gene_symbol"]``, then ``var_names``; missing or blank symbols fail and case
is preserved. Repeated symbols are consolidated by summing their raw counts. Each
cell is then scaled to 10,000 counts across all input genes and log1p-transformed
*before* the official package aligns genes to the model's 5,055-gene panel. Cells
with zero total counts stay all-zero. Cell barcodes (``obs_names``) must be
nonblank and unique, and are preserved. Zero overlap with the model gene panel
fails; low overlap only triggers the package's own warning.

Build and run
-------------

Build on an amd64 Docker host, from the project root::

    docker build --platform linux/amd64 \
        -t run-pan-human-azimuth:local \
        cellarium/cas/benchmarking/comparators/pan_human_azimuth

    docker run --rm --platform linux/amd64 \
        -v "$CACHE_DIR:/home/appuser/.cache/panhumanpy" \
        -v "$DATA_DIR:/data" \
        run-pan-human-azimuth:local \
        --input-h5ad /data/sample.h5ad \
        --output-dir-template '/data/pan_human_azimuth_{stem}' \
        --batch-size 256

Repeat ``--input-h5ad`` to annotate several files in one invocation; model objects
are reloaded for each input.

Options:

* ``--input-h5ad``: required, repeatable, existing file.
* ``--output-dir-template``: required, must contain ``{stem}`` (the input file name
  without extension).
* ``--batch-size``: positive integer, default ``8192``; cells per inference minibatch.

Each output directory must be new or empty. Existing nonempty output directories,
and repeated inputs that resolve to the same directory, fail before inference
without overwriting anything. Missing inputs and invalid options also fail before
model weights are loaded.

Outputs
-------

Each accepted input gets one directory containing:

``pan_human_azimuth_predictions.csv``
    Index is the original cell barcodes, in input order (no index name). Columns,
    in order: ``azimuth_fine``, ``azimuth_medium``, ``azimuth_broad``,
    ``full_hierarchical_labels``, ``final_level_labels``,
    ``final_level_confidence``, ``full_consistent_hierarchy``. Values are the
    package's native output, unmodified.

``pan_human_azimuth_stats.json``
    Keys: ``input_h5ad``, ``n_cells``, ``n_input_genes`` (source gene columns before
    consolidation), ``n_model_genes``, ``n_matched_genes`` (consolidated input
    symbols found in the model panel), ``overlap_ratio``, ``api_version``
    (installed ``panhumanpy``), ``model_version`` (``v1``).

``azimuth_fine``, ``azimuth_medium``, and ``azimuth_broad`` are three resolutions of
one refined annotation, not ranked alternatives. The refined levels have no native
confidence. ``final_level_confidence`` belongs only to ``final_level_labels``, which
can differ from ``azimuth_fine``. No scores or ranks are synthesized.

Unassigned cells
----------------

Model ``v1`` can label poor-quality profiles (for example low-UMI or ambient-RNA
cells) ``Unassigned``. If any cell has ``Unassigned`` at any of the three levels,
the whole input dataset is rejected and none of its outputs are written::

    Pan-Human Azimuth produced Unassigned annotations for 3 cells in '/data/sample.h5ad' (first 3 barcodes): ['AAAC', 'AAAG', 'AACT']

Barcodes are counted once per cell and at most five are shown. A multi-input run
stops at the first failed dataset; outputs of datasets completed earlier are kept.

Python use
----------

::

    from cellarium.cas.benchmarking.comparators.pan_human_azimuth.run import run_pan_human_azimuth

    stats = run_pan_human_azimuth(
        input_h5ad=("/data/sample.h5ad",),
        output_dir_template="/data/pan_human_azimuth_{stem}",
        batch_size=256,
    )

This needs the same runtime dependencies as the image (``panhumanpy==1.0.0``,
``click``, TensorFlow, and so on) and reloads model objects for each input.

Crosswalk
---------

Supply a Pan-Human Azimuth-specific CSV with a label column, ``cl_id``, and
optionally ``cl_label``. This repository does not ship one. The mapped vocabulary
is the union of the ``azimuth_fine``, ``azimuth_medium``, and ``azimuth_broad``
labels. Require every label at every level to appear in the crosswalk, ``cl_id``
values to exist in the benchmark ``ontology_resource.json``, and never
fuzzy-match labels.

Convert and benchmark
---------------------

Run conversion from the project root::

    from cellarium.cas.benchmarking.comparators.pan_human_azimuth.helpers.pan_human_azimuth_to_cas_annotation import (
        pan_human_azimuth_to_cas_annotation,
    )

    pan_human_azimuth_to_cas_annotation(
        predictions_csv_path="/data/pan_human_azimuth_sample/pan_human_azimuth_predictions.csv",
        h5ad_path="/data/sample.h5ad",
        crosswalk_csv_path="/data/pan_human_azimuth_crosswalk.csv",
        crosswalk_label_col="pan_human_label",
        crosswalk_cl_id_col="cl_id",
        crosswalk_cl_label_col="cl_label",
        ontology_resource_path="/data/ontology_resource.json",
        output_dir="/data/pan_human_azimuth_annotate_dir",
    )

The helper writes ``inferred_labels.csv`` with three ranks, as in the existing
Azimuth comparator: rank 1 = ``azimuth_fine``, rank 2 = ``azimuth_medium``,
rank 3 = ``azimuth_broad`` (label, name, and score columns, rows in H5AD order).
Ranks are granularity levels, not alternative guesses; a parent term can repeat
a child's term and is kept. Refined levels have no native confidence, so
``cas_cell_type_score_k`` is empty, not invented. It also writes ``metadata.json``
(model group ``pan_human_azimuth_v1``), ``ontology_response.json``, and
``ontology_resource.json``. It rejects missing barcodes, crosswalk labels, and
ontology terms. Benchmark with ``--inferred-label cas_cell_type_name_1`` and
``--f-measure-top-k 3``.
