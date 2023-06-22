# Cell Annotation Service Client
Python Client libraries for Cell Annotation Service

# Installation
```
$ pip install git+https://github.com/broadinstitute/cell-annotation-service-client.git@fg-data-reading-validation
```
# Usage
To use CAS services create a client instance with an API token:

```python3
from cas_cli import CASClient

api_token = "a_very_long_string_with_some_symbols"
cas = CASClient(api_token=api_token)
```

## Annotation
You can annotate 10x h5 matrices from local disk:
```python3

response = cas.annotate_10x_h5_file(filepath="your_path_to_local_h5_file.h5")
```
or anndata file from local disk:
```python3
response = cas.annotate_anndata_file(filepath="your_path_to_local_h5_file.h5ad")
```
or a previously loaded (unnormalized) anndata object:
```python3
import anndata


adata = anndata.read("you_anndata_file.h5ad")
response = cas.annotate_anndata(adata)
```
