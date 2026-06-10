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

The HiClass literature calls these sets ``alpha_i`` and ``beta_i``. In code and
notebook discussion, it is usually clearer to name them by meaning:

* ``predicted_term_sets[i]`` is the set of unique non-empty
  ``cell_type_ontology_term_id`` values in ``response.data[i].matches``; this is
  ``alpha_i`` in the formula.
* ``ground_truth_ancestor_sets[i]`` is the ontology ancestor set for
  ``ground_truths[i]``, including the ground-truth class itself; this is ``beta_i``
  in the formula.

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


Aggregating benchmark outputs by tissue, donor, assay, or another group
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

If you need group-level metrics later, save class-level outputs and aggregate from
those rows. Do not average per-sample summary rows unless you explicitly want every
sample to have equal weight regardless of cell count or class composition.

For flat and hierarchical F-measure, the robust aggregation pattern is:

1. Join class-level benchmark rows to sample metadata.
2. Pool ``tp``, ``fp``, ``fn``, and ``support`` per group and ground-truth class.
3. Recompute per-class precision, recall, and F1.
4. Derive micro, macro, and macro-weighted summaries from the pooled class rows.

Use this when your samples differ by tissue, donor, assay, perturbation, or batch.
The same notebook pattern works for both flat and hierarchical outputs; only
``metric_prefix`` changes.

Prepare grouped class-level rows
--------------------------------

.. code-block:: python

    import pandas as pd

    # Class-level output from either:
    # - flat: cellarium-cas benchmark flat --save-class-level
    # - hierarchical: cellarium-cas benchmark hierarchical-f-measure --save-class-level
    class_level_df = pd.read_csv("flat_benchmark_class_level.csv")
    # class_level_df = pd.read_csv("hierarchical_f_measure_class_level.csv")

    # One row per sample / annotate output directory.
    # Required key: annotate_dir. Other columns are arbitrary grouping variables.
    sample_metadata = pd.read_csv("sample_metadata.csv")
    # columns might include: annotate_dir, tissue, donor, assay

    df = class_level_df.merge(sample_metadata, on="annotate_dir", how="left")

    # Pick the benchmark family and grouping columns.
    metric_prefix = "flat"          # or "hierarchical"
    group_cols = ["tissue", "assay"]  # examples: ["tissue"], ["donor"], ["tissue", "donor"]

    pooled = (
        df.groupby(group_cols + ["ground_truth_class"], as_index=False)
        .agg(
            tp=("tp", "sum"),
            fp=("fp", "sum"),
            fn=("fn", "sum"),
            support=("support", "sum"),
        )
    )

    precision_col = f"{metric_prefix}_precision"
    recall_col = f"{metric_prefix}_recall"
    f1_col = f"{metric_prefix}_f1"

    pooled[precision_col] = pooled["tp"] / (pooled["tp"] + pooled["fp"])
    pooled[recall_col] = pooled["tp"] / (pooled["tp"] + pooled["fn"])
    pooled[f1_col] = (
        2 * pooled[precision_col] * pooled[recall_col]
        / (pooled[precision_col] + pooled[recall_col])
    )
    pooled = pooled.fillna(0.0)

``pooled`` is the grouped class-level table. It is often the most useful table for
inspection because it shows which classes drive each group's performance.

Compute group-level micro metrics
---------------------------------

Micro metrics pool all counts inside a group before division. They answer: “How well
does the model do for a random cell/node-count contribution in this group?” Abundant
classes contribute more.

.. code-block:: python

    micro = (
        pooled.groupby(group_cols, as_index=False)
        .agg(
            tp=("tp", "sum"),
            fp=("fp", "sum"),
            fn=("fn", "sum"),
            n_cells=("support", "sum"),
        )
    )

    micro[f"micro_{metric_prefix}_precision"] = micro["tp"] / (micro["tp"] + micro["fp"])
    micro[f"micro_{metric_prefix}_recall"] = micro["tp"] / (micro["tp"] + micro["fn"])
    micro[f"micro_{metric_prefix}_f1"] = (
        2 * micro[f"micro_{metric_prefix}_precision"] * micro[f"micro_{metric_prefix}_recall"]
        / (micro[f"micro_{metric_prefix}_precision"] + micro[f"micro_{metric_prefix}_recall"])
    )
    micro = micro.fillna(0.0)

Compute group-level macro metrics
---------------------------------

Macro metrics average grouped class metrics equally. They answer: “How well does the
model do across cell types in this group if every class matters equally?” Rare classes
have the same influence as abundant classes.

.. code-block:: python

    macro = (
        pooled.groupby(group_cols, as_index=False)
        .agg(
            **{
                f"macro_{metric_prefix}_precision": (precision_col, "mean"),
                f"macro_{metric_prefix}_recall": (recall_col, "mean"),
                f"macro_{metric_prefix}_f1": (f1_col, "mean"),
            }
        )
    )

Compute group-level macro-weighted metrics
------------------------------------------

Macro-weighted metrics still compute scores per class first, then weight each class by
its support inside the group. They answer: “How well does the model do across classes,
weighted by how common each class is in this group?”

.. code-block:: python

    pooled["weight"] = pooled["support"] / pooled.groupby(group_cols)["support"].transform("sum")

    weighted = (
        pooled.assign(
            weighted_precision=pooled["weight"] * pooled[precision_col],
            weighted_recall=pooled["weight"] * pooled[recall_col],
            weighted_f1=pooled["weight"] * pooled[f1_col],
        )
        .groupby(group_cols, as_index=False)
        .agg(
            **{
                f"macro_weighted_{metric_prefix}_precision": ("weighted_precision", "sum"),
                f"macro_weighted_{metric_prefix}_recall": ("weighted_recall", "sum"),
                f"macro_weighted_{metric_prefix}_f1": ("weighted_f1", "sum"),
            }
        )
    )

Combine the group summaries
---------------------------

.. code-block:: python

    group_summary = micro.merge(macro, on=group_cols).merge(weighted, on=group_cols)
    group_summary.to_csv(f"{metric_prefix}_benchmark_by_group.csv", index=False)
    pooled.to_csv(f"{metric_prefix}_benchmark_by_group_class.csv", index=False)

Interpretation checklist:

* Use **micro** when cell abundance should dominate.
* Use **macro** when every ground-truth class should matter equally.
* Use **macro-weighted** when you want class-aware metrics but still want common
  classes to contribute proportionally to their support.
* For hierarchical F-measure, ``tp``, ``fp``, and ``fn`` are hierarchical ontology
  set-count totals. For flat metrics, they are pooled flat one-vs-rest counts. The
  aggregation pattern is the same, but the meaning of the counts differs.

Saving results
++++++++++++++

Both summary DataFrames have the same shape (one row), so they can be concatenated and
saved together:

.. code-block:: python

    import pandas as pd

    results = pd.concat([summary_ontology, summary_flat], axis=1)
    results.to_csv("benchmark_results.csv", index=False)
