# Cell Annotation Service Client
Python Client libraries for Cell Annotation Service

# Installation
```
$ pip install git+https://github.com/broadinstitute/cell-annotation-service-client.git@fg-annotate
```
# Usage
To use CAS services create a client instance:
```
from casp_cli import service


cli = service.CASPClientService()
```
CAS cli takes anndata.AnnData as an input to proceed with further analysis

## Annotation
```
import anndata


adata = anndata.read("you_anndata_file.h5ad")
response = cli.annotate_anndata(adata)
```
