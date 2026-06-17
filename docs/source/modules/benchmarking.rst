:tocdepth: 3

Benchmarking
############

The ``cellarium-cas`` package ships benchmarking modules that let you evaluate
annotation quality against a labelled reference dataset.  Two metric families are
supported:

* **F-measure** — standard precision, recall, and F1 computed from a flat confusion
  matrix.  Micro metrics treat the problem globally; macro metrics average F1 across
  classes with nonzero support.
* **Hierarchical F-measure** — hierarchical precision, recall, and F1 following the
  Kiritchenko et al. approach.  Predictions are evaluated using ontology ancestor sets
  rather than exact label matches, so partial credit is awarded for predictions that
  are close in the Cell Ontology tree.

The CLI writes one ``cm_raw_k{n}/`` directory per k value (1 … ``--f-measure-top-k``).
Hierarchical F-measure always uses ``cm_raw_k1/`` (top-1 predictions).  Flat F-measure
computes metrics for every k and stores them in a single wide-format CSV with columns
suffixed ``_k{n}`` (e.g. ``f1_micro_k1``, ``f1_micro_k2``).

Prerequisites
+++++++++++++

Install the benchmarking extras::

    pip install "cellarium-cas[benchmark]"

Preparing annotate output directories
++++++++++++++++++++++++++++++++++++++

Run ``cellarium-cas annotate`` with ``--infer-labels`` and ``--save-ontology-resource``
to produce the files required for benchmarking::

    cellarium-cas annotate \
        --input-path labelled_dataset.h5ad \
        --output-dir ./run_model_a_sample_1 \
        --cas-api-token $CAS_API_TOKEN \
        --infer-labels \
        --save-ontology-resource

Each annotate output directory must contain:

* ``inferred_labels.csv`` — predicted cell type labels (``cas_cell_type_name_1``
  holds the most granular CL term ID).
* ``ontology_resource.json`` — the Cell Ontology resource used during annotation.
* ``metadata.json`` — provenance (``input_path``, ``model_name``, …).
* ``ontology_response.json`` — the raw annotation response.

The ground-truth labels are read from the original ``.h5ad`` at
``metadata["input_path"]``.  That column must contain CL ontology term IDs
(e.g. ``"CL:0000540"``) for hierarchical metrics.
Use ``cas_cell_type_name_k`` when the ground-truth column contains CL term IDs;
use ``cas_cell_type_label_k`` only if the ground-truth column and benchmark label
universe are human-readable labels too.


Running the full pipeline via CLI
++++++++++++++++++++++++++++++++++

The ``cellarium-cas benchmark all`` command runs all four pipeline steps in one
call::

    cellarium-cas benchmark all \
        --annotate-dirs ./annotate_outputs \
        --output-dir    ./benchmark_results \
        --gt-label      cell_type_ontology_term_id \
        --inferred-label cas_cell_type_name_1 \
        --f-measure-top-k 3

``--annotate-dirs`` may be a directory whose immediate subdirectories are annotate
output dirs, or a ``.txt`` file listing one path per line.

``--f-measure-top-k`` computes flat F-measure for every k from 1 to the given value.
All k results appear as separate column sets (suffixed ``_k{n}``) in a single wide-format
CSV — no separate output files per k.


After the command completes, ``benchmark_results/`` contains::

    cm_raw_k1/                       # per-sample sparse confusion matrices (top-1 predictions)
    cm_raw_k2/                       # per-sample sparse confusion matrices (top-2 hit predictions)
    ...                              # one directory per k up to --f-measure-top-k
    cm_aggregate_k1/                 # per-model aggregated confusion matrices (top-1)
    ...                              # one directory per k
    f_measure_per_sample.csv         # F-measure per annotate run (wide format; columns suffixed _k{n})
    f_measure_per_group.csv          # F-measure per model (aggregated; wide format)
    hierarchical_f_measure_per_sample.csv
    hierarchical_f_measure_per_group.csv

Running steps individually
+++++++++++++++++++++++++++

The pipeline can also be run step by step, which is useful when you want to
re-compute only the metric CSVs after adding new runs::

    # Step 1 — build per-sample confusion matrices
    cellarium-cas benchmark confusion-matrix \
        --annotate-dirs ./annotate_outputs \
        --output-dir    ./benchmark_results \
        --gt-label      cell_type_ontology_term_id \
        --inferred-label cas_cell_type_name_1 \
        --f-measure-top-k 3

    # Step 2 — aggregate by model name
    cellarium-cas benchmark aggregate \
        --output-dir ./benchmark_results

    # Step 3 — F-measure CSVs
    cellarium-cas benchmark f-measure \
        --output-dir ./benchmark_results

    # Step 4 — hierarchical F-measure CSVs
    cellarium-cas benchmark hierarchical \
        --output-dir ./benchmark_results

Output columns
++++++++++++++

**f_measure_per_sample.csv** and **f_measure_per_group.csv**:

Metric columns are repeated for each k value and suffixed with ``_k{n}`` (e.g.
``tp_k1``, ``f1_micro_k2``).  Identity columns have no suffix.

.. list-table::
   :header-rows: 1

   * - Column
     - Description
   * - ``model_name`` / ``group_name``
     - Model name (per-sample has ``model_name`` + ``test_sample``; per-group has ``group_name``).
   * - ``tp_k{n}``
     - Global true positives (trace of the confusion matrix) for top-k hit predictions.
   * - ``fp_k{n}``
     - Global false positives (total − trace).
   * - ``fn_k{n}``
     - Global false negatives (total − trace).
   * - ``precision_micro_k{n}``
     - Micro precision (= accuracy for single-label multiclass).
   * - ``recall_micro_k{n}``
     - Micro recall.
   * - ``f1_micro_k{n}``
     - Micro F1.
   * - ``f1_macro_k{n}``
     - Macro F1 averaged over classes with nonzero support.
   * - ``precision_macro_k{n}``
     - Macro precision averaged over classes with nonzero support.
   * - ``recall_macro_k{n}``
     - Macro recall averaged over classes with nonzero support.
   * - ``precision_weighted_k{n}``
     - Support-weighted precision (weighted by per-class true support).
   * - ``recall_weighted_k{n}``
     - Support-weighted recall.
   * - ``f1_weighted_k{n}``
     - Support-weighted F1.

**hierarchical_f_measure_per_sample.csv** and **hierarchical_f_measure_per_group.csv**:

.. list-table::
   :header-rows: 1

   * - Column
     - Description
   * - ``h_tp``
     - Global hierarchical true positives.
   * - ``h_fp``
     - Global hierarchical false positives.
   * - ``h_fn``
     - Global hierarchical false negatives.
   * - ``h_precision_micro``
     - Micro hierarchical precision.
   * - ``h_recall_micro``
     - Micro hierarchical recall.
   * - ``h_f1_micro``
     - Micro hierarchical F1.
   * - ``h_f1_macro``
     - Macro hierarchical F1 (per-class, averaged over true classes with nonzero support).
   * - ``h_precision_macro``
     - Macro hierarchical precision (per-node, averaged over nodes with nonzero true support).
   * - ``h_recall_macro``
     - Macro hierarchical recall (per-node, averaged over nodes with nonzero true support).
   * - ``h_precision_weighted``
     - Support-weighted hierarchical precision (weighted by per-node hierarchical true support).
   * - ``h_recall_weighted``
     - Support-weighted hierarchical recall.
   * - ``h_f1_weighted``
     - Support-weighted hierarchical F1.

Using the Python API directly
++++++++++++++++++++++++++++++

Low-level functions operate on arrays and confusion matrices::

    import scipy.sparse
    from cellarium.cas.benchmarking import (
        build_confusion_matrix,
        compute_f_measure_from_cm,
        compute_hierarchical_f_measure_from_cm,
    )
    from cellarium.cas.postprocessing.cell_ontology.cell_ontology_cache import CellOntologyCache

    # Load ontology cache (provides ancestor mapping)
    import json
    with open("ontology_resource.json") as f:
        resource = json.load(f)
    cache = CellOntologyCache(resource)

    # Build a confusion matrix from label lists
    label_order = resource["cl_names"]
    cm = build_confusion_matrix(y_true, y_pred, label_order)

    # Standard F-measure
    f_metrics = compute_f_measure_from_cm(cm)

    # Hierarchical F-measure
    h_metrics = compute_hierarchical_f_measure_from_cm(
        cm, label_order, cache
    )

The high-level pipeline functions mirror the CLI steps::

    from cellarium.cas.cli._benchmark_impl import (
        run_confusion_matrix_step,
        run_aggregate_step,
        run_f_measure_step,
        run_hierarchical_f_measure_step,
    )

    run_confusion_matrix_step("./annotate_outputs", "./bench", "cell_type_ontology_term_id", "cas_cell_type_name_1")
    run_aggregate_step("./bench")
    run_f_measure_step("./bench")
    run_hierarchical_f_measure_step("./bench")

Aggregating custom groups
+++++++++++++++++++++++++

The CLI ``aggregate`` command groups confusion matrices by ``model_name`` automatically.
For custom groupings (e.g. by assay type or tissue), use
``aggregate_confusion_matrices`` directly
in a notebook::

    from cellarium.cas.benchmarking.confusion_matrix import (
        load_confusion_matrix,
        aggregate_confusion_matrices,
        save_confusion_matrix,
    )

    cms = []
    for run_dir in my_custom_group_dirs:
        cm, _ = load_confusion_matrix(run_dir)
        cms.append(cm)

    agg = aggregate_confusion_matrices(cms)

Azimuth integration
+++++++++++++++++++

The ``cellarium.cas.benchmarking.azimuth`` helpers convert Azimuth annotation outputs
into CAS-compatible annotate directories so they can be evaluated with the same pipeline.
See ``cellarium.cas.benchmarking.azimuth.helpers.azimuth_to_cas_annotation`` for usage.

**Important:** when ``level_specs`` is auto-detected, Azimuth levels are ordered
**most granular first** (rank 1 = finest level), matching CAS convention.  If you
pass explicit ``level_specs``, list them most-granular first.
