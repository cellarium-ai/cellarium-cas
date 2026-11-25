:tocdepth: 3

Usage
#####

API Token and the CAS Public Beta Program
-----------------------------------------
In order to use the CAS API, you will need to join the CAS Public Beta program. To obtain your unique API token to join
the public beta program, please navigate to the CAS webpage at `cellarium.ai <https://cellarium.ai/tool/cellarium-cell-annotation-service-cas/>`_,
scroll to the bottom of the page, and `sign up <https://cellarium.ai/cell-annotation-service-cas-access/>`_. We will contact
you with your unique API key as soon as the public beta is available.

Initialization
--------------
Create a **CASClient** instance with your API token::

    from cellarium.cas import CASClient

    api_token = "api_token_string_goes_here"
    cas = CASClient(api_token=api_token)

Annotation
----------

10x Cell Ranger h5 matrices
+++++++++++++++++++++++++++
You can annotate 10x Cell Ranger h5 matrices from your local disk::

    response = cas.annotate_10x_h5_file(filepath="your_path_to_local_h5_file.h5")

Anndata File from the disk
++++++++++++++++++++++++++
::

    response = cas.annotate_anndata_file(filepath="your_path_to_local_h5_file.h5ad")

Anndata File instance
+++++++++++++++++++++
::

    import anndata


    adata = anndata.read("you_anndata_file.h5ad")
    response = cas.annotate_anndata(adata)

Ontology-aware annotation
--------------------------

Use the ontology-aware strategy to get per-cell relevance scores across the Cell Ontology hierarchy::

    cas_ontology_aware_response = cas.annotate_matrix_cell_type_ontology_aware_strategy(
        matrix=adata,
        chunk_size=500,
        feature_ids_column_name="gene_ids",
        feature_names_column_name="index",
    )

Insert the response into the AnnData object. This call fetches the required cell ontology
resource from the backend and stores scores in ``adata.obsm['cas_cl_scores']`` and metadata
in ``adata.uns['cas_metadata']``::

    cas.insert_ontology_aware_response(cas_ontology_aware_response, adata)

Best cell type label assignment
---------------------------------

The cell type call is selected with a curated L1-L4 granularity lookup. The ``CASClient`` fetches the granularity resource from the backend using the same
ontology resource name as the cell ontology resource, and caches the resulting
``CellGranularityCache`` on the client.

Assign the most granular cell type call to individual cells::

    cas.compute_annotation_levels_single(
        adata=adata,
        min_acceptable_score=0.2,
        obs_prefix="cas_annotation",
    )

Or assign calls per cluster::

    cas.compute_annotation_levels_cluster(
        adata=adata,
        min_acceptable_score=0.2,
        cluster_label_obs_column="cluster_label",
        obs_prefix="cas_annotation_cluster",
    )

Or smooth scores over each cell's local kNN neighborhood while keeping one result per cell::

    cas.compute_annotation_levels_knn(
        adata=adata,
        min_acceptable_score=0.2,
        k_neighbors=30,
        obs_prefix="cas_annotation_knn",
    )

This reuses ``adata.obsp['distances']`` when a kNN graph is already present. Otherwise,
it computes one with ``scanpy.pp.neighbors`` by default. Set
``compute_neighbors_if_missing=False`` to require a precomputed graph. The graph must
contain at least ``k_neighbors`` non-self neighbors per cell.

The lower-level postprocessing functions still accept an externally constructed
``CellGranularityCache``. The client fetches and caches the resource automatically; for
offline testing, preload ``client._granularity_cache[resource_name]`` before calling the
client methods.

Candidates below ``min_acceptable_score`` or without a curated granularity mapping are
excluded. The single selected term maximizes curated depth, then Cell Ontology longest-path depth,
then evidence score, with CL term ID as a deterministic tie-breaker. Missing ontology path depths
count as zero.

The selected term's curated mapping is written as four granularity positions, from most to least
granular (position ``1`` is the lookup's ``l4``, position ``4`` its ``l1``):

- ``{obs_prefix}_label_{n}``: readable name of the term at position ``n``.
- ``{obs_prefix}_name_{n}``: its CL ID.
- ``{obs_prefix}_score_{n}``: its aggregated evidence score.

By default (``fill_missing_levels=True``), missing L4 repeats L3 only when L3 exists.
Missing L3 never inherits L2; missing L4 with no L3 stays null (N/A).
Pass ``fill_missing_levels=False`` to leave L4 null too. When no eligible term
exists, all columns are ``None`` rather than a fabricated root-cell prediction.

Deprecated top-k API
~~~~~~~~~~~~~~~~~~~~

``compute_most_granular_top_k_calls_single``, ``compute_most_granular_top_k_calls_cluster`` and
``compute_most_granular_top_k_calls_knn`` emit a ``DeprecationWarning``. The superseded
``get_most_granular_top_k_calls`` selector remains available without a warning.
Legacy methods keep their original behavior: the top ``top_k`` terms ranked by Cell Ontology
longest (or, with ``use_shortest_path=True``, shortest) path depth then score, written as
``{obs_prefix}_{score,name,label}_{k}`` where ``k`` is the candidate rank, with a root-cell prediction
padding missing ranks. Their default prefixes are ``cas_cell_type``, ``cas_cell_type`` and
``cas_knn_cell_type``, which differ from the ``compute_annotation_levels_*`` defaults, so both can
be run side by side. Use ``compute_annotation_levels_*`` instead.

Interactive visualization
--------------------------

Launch the interactive circular tree / UMAP Dash application::

    from cellarium.cas.visualization import CASCircularTreePlotUMAPDashApp

    CASCircularTreePlotUMAPDashApp(
        adata=adata,
        cas_client=cas,
        cluster_label_obs_column="cluster_label",
    ).run(port=8050)