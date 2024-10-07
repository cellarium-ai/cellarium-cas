:tocdepth: 3


Cellarium Cell Annotation Service
#################################

**What is Cell Annotation Service (CAS)?** CAS is a cloud-native software for rapid querying of single-cell omics data leveraging information-rich and compact vector representations and approximate nearest neighbor search engienes. This documentation refers to the client library for CAS, which provides a simple interface for querying the CAS backend and visualizing the results.

**Hoes does CAS work?** CAS leverages a vast repository of public single-cell data, enriched with metadata such as cell type, disease, tissue, etc. to provide robust annotations for newly generated datasets. As more reference data becomes available and added to the Cellarium DataStore, CAS aims to contextualize new datasets against this accumulated knowledge, addressing a critical bottleneck in single-cell analysis. This approach is akin to how search engines have indexed and made vast amounts of internet data findable and easily accessible. Understanding the mechanisms underlying diseases critically depends on identifying similarities and differences between cellular measurements in various contexts. CAS simplifies this process by identifying similar reference cells and then providing streamlined, community-consensus annotations directly based the context of each queried cell.

.. toctree::
   :maxdepth: 1
   :caption: General Info

   modules/installation
   modules/usage
   modules/workflow
   modules/changelog

.. toctree::
   :maxdepth: 1
   :caption: Codebase Documentation

   automodules/client
   automodules/visualization


Related Projects
================
    * `Cellarium ML Library <https://cellarium-ai.github.io/cellarium-ml/index.html>`_
    * `Cellarium Cloud (Backend) <https://github.com/cellarium-ai/cellarium-cloud>`_

