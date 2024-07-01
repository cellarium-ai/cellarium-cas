:tocdepth: 3

Workflow
########

Token Validation
++++++++++++++++

When you initialize your **CASClient**, it will first check the validity of your API token with the CAS backend.
If your token is valid and there are no other issues, it will print a list of available models for you to use for
annotating your cell data.

Data Validation
+++++++++++++++

Next, when you pass your data to one of the client object's annotate methods, it will validate the format of your
data to ensure it can be properly processed by CAS. At the moment, this means checking if the matrix corresponds to
a CAS schema before sending it to the backend.

Data Sanitization
+++++++++++++++++

If there are any issues during validation that can be resolved with sanitization, the client will do one of the following:

    1. Fill missing features with zeros.
    2. Omit extra features.

Check your quota
++++++++++++++++

Your CAS account has a quota for the number of cells that you can annotate each week.  If the data you have provided
would exceed that quota, the client will raise an exception and inform you of the issue.

The default quota value for a new account is 50,000 cells per week.  You can increase your quota by filling out the
feedback form that is linked in the output when you complete a call to annotate your data.

You can check your total and remaining quota, along with its reset date, by calling the following method on **CASClient**::

    cas.print_user_quota()

Data Transfer to the Backend
++++++++++++++++++++++++++++
If all those checks pass without issue, the **CASClient** will split your input data matrices into smaller shards and
asynchronously send them to the backend to be annotated.  The processing is done in shards as opposed to one large
batch to maximize processing speed.
