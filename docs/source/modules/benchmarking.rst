:tocdepth: 3

Benchmarking
############

The ``cellarium-cas`` package ships complementary benchmarking modules that let you
evaluate annotation quality against a labelled reference dataset:

* **Ontology-aware benchmarking** (``cellarium.cas.benchmarking.compute_ontology_aware_metrics``) —
  measures how close the predicted Cell Ontology (CL) terms are to the ground truth within the
  ontology graph, using a hop-neighbourhood metric.
* **Hierarchical F-measure benchmarking** (``cellarium.cas.benchmarking.compute_hierarchical_f_measure_metrics``) —
  computes HiClass-compatible micro hP/hR/hF and a support-weighted per-ground-truth-class hP/hR/hF.
* **Flat benchmarking** (``cellarium.cas.benchmarking.compute_flat_metrics``) —
  taxonomy-agnostic top-k flat precision/recall/F1 with the same summary/class-level shape as hierarchical F-measure.

Prerequisites
+++++++++++++

Install the benchmarking extras:

.. code-block:: bash

    pip install "cellarium-cas[benchmark]"

Step 1 — Initialise the client and load your data
++++++++++++++++++++++++++++++++++++++++++++++++++

Your input dataset must have ground-truth cell type annotations stored in ``adata.obs``.
For ontology-aware benchmarking the column must contain `Cell Ontology term IDs
<https://www.ebi.ac.uk/ols/ontologies/cl>`_ (e.g. ``"CL:0000540"``).
For flat benchmarking it must contain human-readable cell type labels that match the
labels CAS will predict (e.g. ``"neuron"``).

.. code-block:: python

    import anndata
    from cellarium.cas import CASClient

    cas = CASClient(api_token="<your-api-token>")

    # Load a labelled dataset — ground truth must be in adata.obs
    adata = anndata.read_h5ad("labelled_dataset.h5ad")

    # Column names used in the examples below
    GT_CL_TERM_COLUMN = "cell_type_ontology_term_id"  # e.g. "CL:0000540"
    GT_LABEL_COLUMN   = "cell_type"                   # e.g. "neuron"

Step 2 — Run CAS annotation
++++++++++++++++++++++++++++

Both benchmarking paths start from the ontology-aware annotation strategy, which returns
per-cell relevance scores across the Cell Ontology hierarchy.

.. code-block:: python

    response = cas.annotate_matrix_cell_type_ontology_aware_strategy(
        matrix=adata,
        chunk_size=500,
        feature_ids_column_name="gene_ids",   # column in adata.var holding Ensembl IDs
        feature_names_column_name="index",    # optional; set to None if not available
    )

``response`` is a :class:`~cellarium.cas.models.CellTypeOntologyAwareResults` object whose
``data`` attribute is aligned row-by-row with ``adata``.

Step 3 — Ontology-aware benchmarking
+++++++++++++++++++++++++++++++++++++

Fetch the ontology resource
---------------------------

The benchmarking function needs the raw ontology graph that CAS used internally.
Retrieve it from the backend with the same model that produced the response:

.. code-block:: python

    resource_name    = cas._resolve_ontology_resource_name(response.model_name)
    ontology_resource = cas.cas_api_service.get_cell_ontology_resource(resource_name)

``ontology_resource`` is a plain dictionary with keys ``cl_names``,
``children_dictionary``, and others — exactly the format expected by
``compute_ontology_aware_metrics``.

Run the metric
--------------

.. code-block:: python

    from cellarium.cas.benchmarking import compute_ontology_aware_metrics

    ground_truths = adata.obs[GT_CL_TERM_COLUMN].tolist()

    summary = compute_ontology_aware_metrics(
        response=response,
        ground_truths=ground_truths,
        ontology_resource=ontology_resource,
        num_hops=4,          # evaluate at hop-0 through hop-4
        cell_level=False,    # True → one row per cell instead of a summary row
    )
    print(summary)

The summary ``DataFrame`` has a single row with columns
``n_cells``, ``hop_0_tpr``, ``hop_0_tnr``, ``hop_0_f1_score``, …
``hop_4_tpr``, ``hop_4_tnr``, ``hop_4_f1_score``.

**Hop semantics**: at hop-N a predicted CL term is considered a *true positive* if it
lies within N bidirectional steps of the ground-truth term in the ontology graph.
Hop-0 is an exact match; each additional hop widens the acceptable neighbourhood by one
parent/child edge.

Cell-level output
-----------------

Pass ``cell_level=True`` to get a row per cell, useful for per-cluster analysis:

.. code-block:: python

    cell_df = compute_ontology_aware_metrics(
        response=response,
        ground_truths=ground_truths,
        ontology_resource=ontology_resource,
        num_hops=4,
        cell_level=True,
    )
    # Columns: query_cell_id, ground_truth, hop_0_tpr, …, hop_4_f1_score

Step 4 — Hierarchical F-measure benchmarking
++++++++++++++++++++++++++++++++++++++++++++

The hierarchical F-measure path treats every cell as a pair of ontology node sets:

* ``alpha_i`` is the set of unique non-empty ``cell_type_ontology_term_id`` values in
  ``response.data[i].matches``.
* ``beta_i`` is the ontology ancestor set for ``ground_truths[i]``, including the
  ground-truth class itself.

.. code-block:: python

    from cellarium.cas.benchmarking import compute_hierarchical_f_measure_metrics

    hierarchical_summary = compute_hierarchical_f_measure_metrics(
        response=response,
        ground_truths=adata.obs[GT_CL_TERM_COLUMN].tolist(),
        ontology_resource=ontology_resource,
    )

The summary ``DataFrame`` has exactly these columns: ``n_cells``,
``micro_hierarchical_precision``, ``micro_hierarchical_recall``,
``micro_hierarchical_f1``, ``macro_hierarchical_precision``,
``macro_hierarchical_recall``, ``macro_hierarchical_f1``,
``macro_weighted_hierarchical_precision``, ``macro_weighted_hierarchical_recall``,
and ``macro_weighted_hierarchical_f1``.

The micro columns use the literature-defined HiClass ``average=\"micro\"`` set-count
ratios. The macro columns are unweighted means of the per-ground-truth-class rows.
The macro-weighted columns use the same class rows, weighted by class support.

When used through ``cellarium-cas benchmark hierarchical-f-measure``, the summary CSV
contains one ``row_type=\"sample\"`` row per annotate output directory and one
``row_type=\"total\"`` row per model. Total rows pool all cells for that model first,
then recompute both micro and macro-weighted metrics; they are not averages of the
sample rows.

Pass ``class_level=True`` to return one row per ground-truth class with columns
``ground_truth_class``, ``support``, ``weight``, ``tp``, ``fp``, ``fn``,
``hierarchical_precision``, ``hierarchical_recall``, and ``hierarchical_f1``. In the
CLI, ``--save-class-level`` writes a single ``hierarchical_f_measure_class_level.csv``
with per-class ``row_type=\"sample\"`` rows and per-model ``row_type=\"total\"`` rows.

Step 5 — Flat (taxonomy-agnostic) benchmarking
++++++++++++++++++++++++++++++++++++++++++++++

The flat benchmarking path compares predicted human-readable labels against ground-truth
labels.  First, insert the response into ``adata`` and derive ranked label predictions.

Insert the response and assign labels
--------------------------------------

.. code-block:: python

    # Writes cas_cl_scores into adata.obsm and metadata into adata.uns
    cas.insert_ontology_aware_response(response, adata)

    # Assigns top-k most granular cell type labels to adata.obs columns:
    # cas_cell_type_label_1, cas_cell_type_label_2, …, cas_cell_type_label_5
    cas.compute_most_granular_top_k_calls_single(
        adata=adata,
        min_acceptable_score=0.1,
        top_k=5,
        obs_prefix="cas_cell_type",
    )

Extract predictions and run the metric
---------------------------------------

.. code-block:: python

    from cellarium.cas.benchmarking import compute_flat_metrics, extract_predictions_from_adata

    # Reads cas_cell_type_label_1 … cas_cell_type_label_5 from adata.obs
    predictions = extract_predictions_from_adata(
        adata,
        column_prefix="cas_cell_type_label",
        top_k=5,
    )

    ground_truths = adata.obs[GT_LABEL_COLUMN].tolist()

    summary = compute_flat_metrics(
        ground_truths=ground_truths,
        predictions=predictions,
        top_k=5,           # evaluate using the effective top-5 prediction
        class_level=False, # True → one row per ground-truth class
    )
    print(summary)

The summary ``DataFrame`` has the same shape as hierarchical F-measure output: one row
with ``n_cells``, ``micro_flat_precision``, ``micro_flat_recall``, ``micro_flat_f1``,
``macro_flat_precision``, ``macro_flat_recall``, ``macro_flat_f1``,
``macro_weighted_flat_precision``, ``macro_weighted_flat_recall``, and
``macro_weighted_flat_f1``. Pass ``class_level=True`` to get the matching class-level
shape: ``ground_truth_class``, ``support``, ``weight``, ``tp``, ``fp``, ``fn``,
``flat_precision``, ``flat_recall``, and ``flat_f1``.

Saving results
++++++++++++++

Both summary DataFrames have the same shape (one row), so they can be concatenated and
saved together:

.. code-block:: python

    import pandas as pd

    results = pd.concat([summary_ontology, summary_flat], axis=1)
    results.to_csv("benchmark_results.csv", index=False)
