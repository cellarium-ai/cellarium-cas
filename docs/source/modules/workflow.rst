:tocdepth: 3

Workflow
#######

Token Validation
++++++++++++++++

**CASClient** first checks the API token validity

Data Validation
+++++++++++++++

Next step is to check the data. At the moment, **CASClient** checks if the matrix corresponds to CAS schema before
sending it to the backend.

Data Sanitization
+++++++++++++++++

Depending on the previous step it would next do one of the following:

    1. Fill missing features with zeros.
    2. Omit extra features.

Data Transfer to Backend
++++++++++++++++++++++++
**CASClient** splits input data matrices into smaller shards and asynchronously send them to backend to achieve the
maximum speed of processing
