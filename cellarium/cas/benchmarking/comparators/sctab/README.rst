scTab Benchmarking
==================

This directory contains the Docker image and helper scripts used to run
`scTab <https://github.com/theislab/scTab>`_ cell type annotations and convert them
into the format expected by ``cellarium-cas benchmark``.  scTab is a cross-tissue
classifier that predicts one of 164 cell types from 19,331 protein-coding gene counts,
normalized to 10,000 counts + log1p.  All evaluation is performed by the same
``cellarium-cas benchmark`` commands used for CAS, ensuring a fair, apples-to-apples
comparison.

scTab is pinned to upstream commit ``5ede7f2``.  The container runs CPU-only inference
with a per-cell top-k of softmax probabilities; the host helpers map those class labels
to Cell Ontology (CL) terms and build the same four output files as ``cellarium-cas
annotate``.

Directory structure::

    Dockerfile                        # builds the scTab Docker image (CPU, python:3.10-slim)
    run.py                            # entrypoint: h5ad → top-k predictions CSV
    helpers/
        map_sctab_to_cas_labels.py    # predictions CSV → inferred_labels.csv + metadata.json
        sctab_to_cas_annotation.py    # predictions CSV → all four annotate files

The shared ``ontology_response.json`` builder lives at
``../shared/cas_labels.py`` and is used by both comparators.


Public Docker Image
-------------------

A pre-built image is planned at::

    us-central1-docker.pkg.dev/dsp-cellarium/dev-resources/run-sctab:0.0.1

The image bundles CPU PyTorch, scTab (commit ``5ede7f2``), and the inference script.
It loads one checkpoint and processes one or more ``.h5ad`` files, writing per-cell
top-k predictions to CSV.

To build the image from source (optional)::

    cd benchmarking/comparators/sctab
    docker build -t run-sctab:local .


Crosswalk Reference
-------------------

scTab class labels are mapped to Cell Ontology (CL) term IDs using the same HRA
crosswalk table as Azimuth::

    https://apps.humanatlas.io/kg-explorer/ctann/azimuth/v1.3

Download the crosswalk CSV from that page.  Open the file and note the exact column
names before running the helpers — typical column names are:

- ``Annotation_Label`` — scTab class label (``crosswalk_label_col``)
- ``CL_ID`` — Cell Ontology term ID (``crosswalk_cl_id_col``)
- ``CL_Label`` — human-readable CL term label (``crosswalk_cl_label_col``)


Step-by-step Usage
------------------

**Step 1 — Pull the Docker image**
::

    docker pull us-central1-docker.pkg.dev/dsp-cellarium/dev-resources/run-sctab:0.0.1

**Step 2 — Run scTab on your .h5ad file**

Mount the directory containing your data and call the image with the checkpoint,
metadata files, and input::

    docker run --rm \
        -v /path/to/data:/data \
        -v /path/to/checkpoints:/checkpoints \
        us-central1-docker.pkg.dev/dsp-cellarium/dev-resources/run-sctab:0.0.1 \
        --model-type tabnet \
        --checkpoint /checkpoints/scTab_tabnet.ckpt \
        --hparams /checkpoints/hparams.yaml \
        --var-parquet /checkpoints/var.parquet \
        --cell-type-parquet /checkpoints/cell_type.parquet \
        --input-h5ad /data/my_dataset.h5ad \
        --output-dir-template /data/sctab_{stem} \
        --top-k 5

``model_type`` is one of ``tabnet``, ``mlp``, ``linear`` and selects the classifier
architecture.  The checkpoint's ``classifier.*`` weights are loaded into that
architecture.  Run the image with ``--help`` for the full flag table (``--chunk-size``,
``--min-gene-overlap``).

Outputs per input (one directory named from ``--output-dir-template`` with ``{stem}``
replaced by the input filename without extension):

- ``sctab_predictions.csv`` — one row per cell (index = barcodes), columns
  ``predicted_label_1, predicted_score_1, …, predicted_label_K, predicted_score_K``
  (softmax probabilities).
- ``sctab_stats.json`` — ``n_cells``, ``n_model_genes``, ``n_matched_genes``,
  ``overlap_ratio``, ``model_type``.

Input gene ids are matched to the model's 19,331-gene vocabulary by version-stripped
Ensembl id; missing model genes are zero-filled.  A run raises if fewer than
``--min-gene-overlap`` (default 0.5) of model genes are present in an input.

**Step 3 — Download the crosswalk**

Download the HRA crosswalk CSV from::

    https://apps.humanatlas.io/kg-explorer/ctann/azimuth/v1.3

Save it locally, e.g. as ``crosswalk.csv``.

**Step 4 — Convert scTab predictions to CAS inferred_labels format**

Run from the project root (``cellarium-cas`` must be installed)::

    from cellarium.cas.benchmarking.comparators.sctab.helpers.map_sctab_to_cas_labels import map_sctab_to_cas_labels

    map_sctab_to_cas_labels(
        predictions_csv_path="/data/sctab_my_dataset/sctab_predictions.csv",
        h5ad_path="/data/my_dataset.h5ad",
        output_dir="/data/sctab_annotate_dir",
        crosswalk_csv_path="crosswalk.csv",
        crosswalk_label_col="Annotation_Label",
        crosswalk_cl_id_col="CL_ID",
        crosswalk_cl_label_col="CL_Label",      # optional but recommended
        ontology_resource_path="/path/to/ontology_resource.json",
        model_name="sctab_tabnet",
        top_k=5,
    )

The helper validates that every predicted label maps through the crosswalk and that
every mapped CL id exists in the ontology resource, then writes ``inferred_labels.csv``
and ``metadata.json`` inside the output directory.  Ranks mapping to the same CL id are
merged by summing their scores.

**Step 5 — Build the CAS-compatible ontology response**

Or run Steps 4 and 5 together with the combined helper (this builds
``ontology_response.json`` and copies ``ontology_resource.json`` using the shared
builder)::

    from cellarium.cas.benchmarking.comparators.sctab.helpers.sctab_to_cas_annotation import sctab_to_cas_annotation

    sctab_to_cas_annotation(
        predictions_csv_path="/data/sctab_my_dataset/sctab_predictions.csv",
        h5ad_path="/data/my_dataset.h5ad",
        crosswalk_csv_path="crosswalk.csv",
        crosswalk_label_col="Annotation_Label",
        crosswalk_cl_id_col="CL_ID",
        ontology_resource_path="/path/to/ontology_resource.json",
        model_name="sctab_tabnet",
        output_dir="/data/sctab_annotate_dir",
        top_k=5,
    )

**Step 6 — Run the benchmark**

After Step 5, the output directory contains all four files required by
``cellarium-cas benchmark`` (``inferred_labels.csv``, ``metadata.json``,
``ontology_response.json``, ``ontology_resource.json``).

Run the full pipeline in one command::

    cellarium-cas benchmark all \
        --annotate-dirs  /data \
        --output-dir     /data/benchmark_results \
        --gt-label       <your_ground_truth_obs_column> \
        --inferred-label cas_cell_type_name_1

Here ``/data`` is the **parent directory** containing ``sctab_annotate_dir`` as a
subdirectory.  Point ``--annotate-dirs`` at the same parent that holds your CAS
annotate output directories to compare all tools in one run.


Notes
-----

- The ``model_name`` field is ``sctab_<model_type>`` — ``sctab_tabnet``,
  ``sctab_mlp``, or ``sctab_linear`` — which appears in the benchmark summary CSVs
  alongside CAS model names.  Each model type forms its own benchmark group.
- ``cellarium-cas`` must be installed in the Python environment used to run the
  helpers (Steps 4 and 5).  The Docker image is only needed for Step 2.
- The helpers require ``anndata``, ``pandas``, and ``scipy`` (already dependencies of
  ``cellarium-cas``).  The image installs torch, lightning, torchmetrics, pyarrow,
  and scTab itself.