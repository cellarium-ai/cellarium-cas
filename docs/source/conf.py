# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import shutil
import glob
import time

from setuptools_git_versioning import get_tag

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Cellarium CAS"
copyright = f"{time.strftime('%Y')}, Cellarium AI Lab"
author = "Cellarium AI Lab"
version = get_tag() or "<no version>"
release = get_tag() or "<no release>"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'nbsphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_substitution_extensions",
    "sphinxcontrib.autodoc_pydantic",
    "IPython.sphinxext.ipython_console_highlighting",
]

exclude_patterns = ['_build', '**.ipynb_checkpoints']

# Provide substitutions for common values
rst_prolog = f"""
.. |project| replace:: {project}
.. |copyright| replace:: {copyright}
.. |author| replace:: {author}
.. |version| replace:: {version}
.. |release| replace:: {release}
"""

rst_epilog = """
.. |br| raw:: html

   <br />
"""

intersphinx_mapping = {
    "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"

nitpicky = True
nitpick_ignore_regex = [
    # Ignore exceptions from nested Pydantic models
    (r'py:.*', r'cellarium\.cas\.models\..*'),
]

# The JSON schema is a bit much in the docs
autodoc_pydantic_model_show_json = False

if not os.path.exists("notebooks"):
    os.makedirs("notebooks")

for src_file in glob.glob("../../notebooks/*.ipynb"):
    shutil.copy(src_file, "notebooks/")
