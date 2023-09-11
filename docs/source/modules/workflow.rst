:tocdepth: 3

Workflow
########

Token Validation
++++++++++++++++

The **CASClient** first checks the validity of the API token.

Data Validation
+++++++++++++++

Next step is to check the data. At the moment, the **CASClient** checks if the matrix corresponds to a CAS schema
before sending it to the backend.

Data Sanitization
+++++++++++++++++

Depending on the previous step it would next do one of the following:

    1. Fill missing features with zeros.
    2. Omit extra features.

Data Transfer to the Backend
++++++++++++++++++++++++++++
The **CASClient** splits input data matrices into smaller shards and asynchronously send them to the backend to achieve
the maximum speed of processing.
