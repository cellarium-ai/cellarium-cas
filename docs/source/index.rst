:tocdepth: 3


Cellarium Cell Annotation Service (CAS) Documentation
#####################################################

Overview
++++++++

**What is the Cell Annotation Service (CAS)?**  
The Cell Annotation Service (CAS) is a cloud-native platform designed for rapid and efficient querying of single-cell omics data. It utilizes compact and information-rich vector representations combined with approximate nearest neighbor search engines to enable seamless exploration of vast single-cell datasets. This documentation focuses on the CAS client library, which offers a user-friendly interface for querying the CAS backend and visualizing the results.

**Hows does CAS work?**
CAS operates by building a vector search index derived from low-dimensional embeddings of a comprehensive repository of publicly available single-cell transcriptomics data. Currently, the CAS reference dataset comprises the entire `CZI CELLxGENE <https://cellxgene.cziscience.com/>`_ data catalog. This reference database includes rich metadata such as cell type, disease state, tissue origin, and more. In future iterations, the reference catalog will be expanded to incorporate additional datasets.

When a user queries CAS with their single-cell transcriptomics data, the system maps the input data into the same vector space as the reference dataset. It then performs an approximate nearest neighbor search to identify similar cells within the reference database. Based on these similar reference cells, CAS generates label summary statistics, providing comprehensive annotations for each query cell. The methodology employed by CAS is comparable to how inverse search engines (such as image search techniques) index and retrieve information from vast amounts of internet data, making it accessible and interpretable. The embeddings used by CAS are generated using distributed machine learning models implemented in the `Cellarium ML Library <https://cellarium-ai.github.io/cellarium-ml/index.html>`_

Understanding the similarities and differences between cellular measurements in various contexts is crucial for unraveling disease mechanisms. CAS streamlines this process by leveraging its reference context to provide community-consensus annotations for the queried cells, helping researchers gain deeper insights into their data. 

Funding
+++++++

Cellarium CAS was co-developed by 10x Genomics and Cellarium AI Lab at the Data Sciences Platform, Broad Institute. The project was funded by 10x Genomics, NIH Grant UM1 MH130966, and additional support from Broad Institute.

Future Plans
++++++++++++

At present, CAS outputs annotations exclusively for cell types. In upcoming updates, we plan to extend the service to include a broader range of informative metadata such as disease status, developmental stage, tissue origin, and other relevant biological contexts, thereby providing users with a more comprehensive annotation framework for their single-cell data.


.. toctree::
   :maxdepth: 1
   :caption: General Usage

   modules/installation
   modules/usage
   modules/workflow
   modules/changelog

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   Quickstart Tutorial (on Google Colab) <https://colab.research.google.com/drive/1m9zgqP5n7E4pGGCg5RjfvlCnS6uqUdSa>
   Quickstart Tutorial (static) <notebooks/quickstart_tutorial.ipynb>

.. toctree::
   :maxdepth: 1
   :caption: Codebase Documentation

   automodules/client
   automodules/visualization


Related Projects
================
    * `Cellarium ML Library <https://cellarium-ai.github.io/cellarium-ml/index.html>`_
    * `Cellarium Cloud (Backend) <https://github.com/cellarium-ai/cellarium-cloud>`_

