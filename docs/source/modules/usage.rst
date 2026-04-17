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

Assign the most granular top-k cell type calls to individual cells::

    cas.compute_most_granular_top_k_calls_single(
        adata=adata,
        min_acceptable_score=0.2,
        top_k=3,
        obs_prefix="cas_cell_type",
    )

Or assign calls per cluster::

    cas.compute_most_granular_top_k_calls_cluster(
        adata=adata,
        min_acceptable_score=0.2,
        cluster_label_obs_column="cluster_label",
        top_k=3,
        obs_prefix="cas_cell_type_cluster",
    )

Interactive visualization
--------------------------

Launch the interactive circular tree / UMAP Dash application::

    from cellarium.cas.visualization import CASCircularTreePlotUMAPDashApp

    CASCircularTreePlotUMAPDashApp(
        adata=adata,
        cas_client=cas,
        cluster_label_obs_column="cluster_label",
    ).run(port=8050)