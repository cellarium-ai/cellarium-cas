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
version = get_tag()
release = get_tag()

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_substitution_extensions",
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

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
