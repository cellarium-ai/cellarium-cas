Azimuth Benchmarking
====================

This directory contains the Docker image and helper scripts used to run
`Azimuth <https://azimuth.hubmapconsortium.org/>`_ cell type annotations and convert
them into the format expected by ``cellarium-cas benchmark``.  All evaluation is
performed by the same ``cellarium-cas benchmark flat`` and
``cellarium-cas benchmark ontology-aware`` commands used for CAS, ensuring a
fair, apples-to-apples comparison.

Directory structure::

    Dockerfile                        # builds the Azimuth Docker image
    run.sh                            # entrypoint: h5ad → Azimuth CSV
    run_azimuth.R                     # R script called by run.sh
    helpers/
        convert_h5ad.py               # converts .h5ad to 10x format for Seurat
        map_azimuth_to_cas_labels.py  # Azimuth CSV → inferred_labels.csv + metadata.json
        build_ontology_response.py    # Azimuth CSV → ontology_response.json


Public Docker Image
-------------------

A pre-built image is publicly available::

    us-central1-docker.pkg.dev/dsp-cellarium/dev-resources/run-azimuth:0.0.1

The image bundles R, Seurat, Azimuth, and a Python conversion script.  It accepts
a ``.h5ad`` file as input, converts it to 10x Genomics sparse format internally,
runs Azimuth, and writes the annotated cell metadata to a CSV file.

To build the image from source (optional)::

    cd benchmarking/azimuth
    docker build -t run-azimuth:local .


Crosswalk Reference
-------------------

Azimuth cell type labels are mapped to Cell Ontology (CL) term IDs using the HRA
crosswalk table for Azimuth v1.3:

    https://apps.humanatlas.io/kg-explorer/ctann/azimuth/v1.3

Download the crosswalk CSV from that page.  Open the file and note the exact column
names before running the helpers — typical column names are:

- ``tool_cell_label`` — Azimuth cell type label (``--crosswalk-azimuth-col``)
- ``cl_id`` — Cell Ontology term ID (``--crosswalk-cl-id-col``)

The crosswalk covers all Azimuth annotation levels in a single table; no separate
file per level is needed.


Step-by-step Usage
------------------

**Step 1 — Pull the Docker image**
::

    docker pull us-central1-docker.pkg.dev/dsp-cellarium/dev-resources/run-azimuth:0.0.1

**Step 2 — Run Azimuth on your .h5ad file**

Mount the directory containing your data and call the image with three positional
arguments: ``<input_h5ad_path> <azimuth_reference> <output_csv_path>``::

    docker run --rm \
        -v /path/to/data:/data \
        us-central1-docker.pkg.dev/dsp-cellarium/dev-resources/run-azimuth:0.0.1 \
        /data/my_dataset.h5ad \
        pbmcref \
        /data/azimuth_output.csv

Replace ``pbmcref`` with the appropriate Azimuth reference (e.g. ``bonemarrowref``,
``lungref``).  See the `Azimuth references page
<https://azimuth.hubmapconsortium.org/references/>`_ for the full list.

The output CSV has one row per cell with columns for each annotation level, for
example ``predicted.celltype.l1``, ``predicted.celltype.l1.score``,
``predicted.celltype.l2``, ``predicted.celltype.l2.score``.

**Step 3 — Download the crosswalk**

Download the HRA crosswalk CSV from::

    https://apps.humanatlas.io/kg-explorer/ctann/azimuth/v1.3

Save it locally, e.g. as ``crosswalk.csv``.

**Step 4 — Convert Azimuth labels to CAS inferred_labels format**
::

    python cellarium/cas/benchmarking/azimuth/helpers/map_azimuth_to_cas_labels.py \
        --azimuth-csv            /data/azimuth_output.csv \
        --h5ad-path              /data/my_dataset.h5ad \
        --output-dir             /data/azimuth_annotate_dir \
        --crosswalk-csv          crosswalk.csv \
        --crosswalk-azimuth-col  tool_cell_label \
        --crosswalk-cl-id-col    cl_id \
        --azimuth-ref-name       pbmcref \
        --level predicted.celltype.l3:predicted.celltype.l3.score \
        --level predicted.celltype.l2:predicted.celltype.l2.score \
        --level predicted.celltype.l1:predicted.celltype.l1.score

List ``--level`` arguments from **most granular to least granular**.  Each value is
``<label_column>:<score_column>`` as they appear in the Azimuth output CSV.

This produces ``inferred_labels.csv`` and ``metadata.json`` inside the output directory.

**Step 5 — Build the CAS-compatible ontology response**

This step requires an ``ontology_resource.json`` file (saved by
``cellarium-cas annotate --save-ontology-resource``) and the
``inferred_labels.csv`` produced in Step 4::

    python cellarium/cas/benchmarking/azimuth/helpers/build_ontology_response.py \
        --inferred-labels-path   /data/azimuth_annotate_dir/inferred_labels.csv \
        --output-path            /data/azimuth_annotate_dir/ontology_response.json \
        --ontology-resource-path /path/to/ontology_resource.json \
        --azimuth-ref-name       pbmcref

**Step 6 — Run the benchmark**

After Steps 4 and 5, ``/data/azimuth_annotate_dir`` contains all three files
required by ``cellarium-cas benchmark``.

Flat benchmark (taxonomy-agnostic)::

    cellarium-cas benchmark flat \
        --annotate-dirs  /data \
        --gt-column-name <your_ground_truth_column> \
        --output-dir     /data/benchmark_results

Ontology-aware benchmark::

    cellarium-cas benchmark ontology-aware \
        --annotate-dirs      /data \
        --gt-cl-column-name  <your_cl_id_column> \
        --output-dir         /data/benchmark_results

Here ``/data`` is the **parent directory** containing ``azimuth_annotate_dir`` as a
subdirectory.  Point ``--annotate-dirs`` at the same parent that holds your CAS
annotate output directories to compare all tools in one run.


Notes
-----

- The ``model_name`` field in outputs is set to ``azimuth_<ref_name>``
  (e.g. ``azimuth_pbmcref``), which appears in the benchmark summary CSVs alongside
  CAS model names.
- ``cellarium-cas`` must be installed in the Python environment used to run the
  helpers (Steps 4 and 5).  The Docker image is only needed for Step 2.
- The helpers require ``anndata``, ``pandas``, and ``scipy`` (already dependencies of
  ``cellarium-cas``).
