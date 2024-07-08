# Cellarium Cell Annotation Service (CAS) Client
This codebase contains the Python client library for using Cellarium Cell Annotation Service (CAS).

# Installation
```
$ pip install cellarium-cas
```
# Usage
To use Cellarium CAS, create a client instance with your API token:

```python3
from cellarium.cas import CASClient

api_token = "a_very_long_string_with_some_symbols"
cas = CASClient(
  api_token=api_token,
  api_url="<optional url to connect to a non-standard CAS server>"
)
```

## Annotation
You can annotate 10x Cell Ranger h5 matrices from local disk:
```python3

response = cas.annotate_10x_h5_file(filepath="your_path_to_local_h5_file.h5")
```
or an anndata file from local disk:
```python3
response = cas.annotate_anndata_file(filepath="your_path_to_local_h5_file.h5ad")
```
or a previously loaded (unnormalized) anndata object:
```python3
import anndata


adata = anndata.read("you_anndata_file.h5ad")
response = cas.annotate_anndata(adata)
```
