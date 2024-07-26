:tocdepth: 3

Usage
#####

API Token and the CAS Early Access Program
------------------------------------------
In order to use the CAS API, you will need to join the Cell Annotation Service Early Access program.  To join the program,
please fill out the form at https://cellarium.ai/cell-annotation-service-cas-early-access/ and we will notify you
when your account is created, and provide you with an API token to use with the CAS API.

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