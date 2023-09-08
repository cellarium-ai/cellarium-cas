:tocdepth: 3

Usage
#####

API Token
---------
You need to request an API token before you can use it.
Please contact our team by mehrtash@broadinstitute.org if you're interested in our tool and we will create an
API token for your needs.

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
You can annotate 10x Cell Ranger h5 matrices from local disk::

    response = cas.annotate_10x_h5_file(filepath="your_path_to_local_h5_file.h5")

Anndata File from disk
++++++++++++++++++++++
::

    response = cas.annotate_anndata_file(filepath="your_path_to_local_h5_file.h5ad")

Anndata File instance
+++++++++++++++++++++
::

    import anndata


    adata = anndata.read("you_anndata_file.h5ad")
    response = cas.annotate_anndata(adata)