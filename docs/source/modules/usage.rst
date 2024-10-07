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