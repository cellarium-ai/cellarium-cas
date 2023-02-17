# Cell Annotation Service Client
Python Client libraries for Cell Annotation Service

# Installation
```
$ pip install git+https://github.com/broadinstitute/cell-annotation-service-client.git
```
# Usage
To use CAS services create a client instance with an API token:
```python3
from casp_cli import service


api_token = "a_very_long_string_with_some_symbols"
cli = service.CASClientService(api_token=api_token)
```
CAS cli takes anndata.AnnData as an input to proceed with further analysis

## Annotation
```
import anndata


adata = anndata.read("you_anndata_file.h5ad")
response = cli.annotate_anndata(adata)
```
