# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import time

from setuptools_git_versioning import get_tag

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Cellarium CAS"
copyright = f"{time.strftime('%Y')}, Cellarium AI"
author = "Cellarium AI"
version = get_tag() or "<no version>"
release = get_tag() or "<no release>"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx_substitution_extensions",
    "sphinxcontrib.autodoc_pydantic",
]

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
