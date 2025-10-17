<img src="https://cellarium.ai/wp-content/uploads/2024/07/cellarium-logo-medium.png" alt="Cellarium Logo" width="180">

# Cellarium Cell Annotation Service (CAS) Client Library
This codebase contains the Python client library for using Cellarium Cell Annotation Service (CAS).

# Installation
The cellarium-cas package officially supports Python versions between 3.7 and 3.12.  We recommend using Python 3.10+.
You can install CAS from PyPI using `pip`. To install the latest version, please run the following command:
```
$ pip install cellarium-cas
```
To install a specific version `version_number`, you can use the following command:
```
$ pip install cellarium-cas==<version_number>
```
If you wish to use visualization features, you can install the package with the visualization extras:
```
$ pip install cellarium-cas[vis]
```

# Obtaining an API Token
You need an API token to use CAS. We are offerring a free public beta program for a limited time to try CAS and explore ways it can enhance your cell biology research. To obtain your unique API token to join the public beta program, please navigate to the CAS webpage at [cellarium.ai](https://cellarium.ai/tool/cellarium-cell-annotation-service-cas/), scroll to the bottom of the page, and [sign up](https://cellarium.ai/cell-annotation-service-cas-access/). We will contact you with your unique API key as soon as the public beta is available.

# Quickstart Tutorial
The fastest way to get started with using CAS is to follow the quickstart tutorial:
[Click here to open the quickstart tutorial on GitHub](notebooks/quickstart_tutorial.ipynb)

It is even easier to go through the quickstart tutorial on Google Colab. Remember, you still need an API key to successfully run through the tutorial:
[Click here to open the quickstart tutorial on Google Colab](https://colab.research.google.com/drive/1m9zgqP5n7E4pGGCg5RjfvlCnS6uqUdSa)

# Documentation
Please visit the project's [ReadTheDocs page](https://cellarium-cas.readthedocs.io/) for additional documentation.

# Citation
If the Cellarium Cell Annotation Service (CAS) contributes to your research, please cite our preprint to help others discover the project. You can reference it as follows:

> Williams, Stephen R., Grab, Fedor, Kamath, Govinda M., Ordabayev, Yerdos, Mellen, Jeff, Roelli, Patrick, Cibulskis, Kristian, Lehnert, Erik, Xie, Fen, Covarrubias, Miguel, Rahman, Nur-Taz, Tickle, Timothy, Erhan, Emre, Malfroy-Camine, Nicolas, Lydon, Kevin, Babadi, Mehrtash, and Delaney, Nigel F. 2025. *Accelerating scRNA-seq Analysis: Automated cell type annotation using representation learning and vector search.* bioRxiv. https://doi.org/10.1101/2025.10.06.680787

BibTeX users can cite CAS with:

```
@article{cellarium_10x_cas,
  author = {Williams, Stephen R. and Grab, Fedor and Kamath, Govinda M. and Ordabayev, Yerdos and Mellen, Jeff and Roelli, Patrick and Cibulskis, Kristian and Lehnert, Erik and Xie, Fen and Covarrubias, Miguel and Rahman, Nur-Taz and Tickle, Timothy and Erhan, Emre and Malfroy-Camine, Nicolas and Lydon, Kevin and Babadi, Mehrtash and Delaney, Nigel F.},
  title = {Accelerating scRNA-seq Analysis: Automated cell type annotation using representation learning and vector search},
  journal = {bioRxiv},
  year = {2025},
  doi = {10.1101/2025.10.06.680787},
  url = {https://www.biorxiv.org/content/10.1101/2025.10.06.680787v1}
}
```
