# Cellarium Cloud Cell Annotation Service (CAS) Client
Python client libraries

# Installation
```
$ pip install git+https://github.com/broadinstitute/cell-annotation-service-client.git
```
# Usage
To use CAS services create a client instance:
```
from cas_cli import service


cli = service.CASClientService()
```
CAS cli takes anndata.AnnData as an input to proceed with further analysis

## Annotation
```
import anndata


adata = anndata.read("you_anndata_file.h5ad")
response = cli.annotate_anndata(adata)
```
