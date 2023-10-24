# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
from setuptools_git_versioning import get_tag

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Cellarium CAS"
copyright = "2023, Cellarium AI"
author = "Cellarium AI"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_rtd_theme",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["**/installation.rst"]

rst_epilog = """
.. |br| raw:: html

   <br />
"""

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]


def substitute_installation_version():
    """
    Generate /modules/installation_.rst from /modules/installation.rst by substituting '|cas_version|' with a tagged
    version of a library.
    This step was necessary because the Sphinx builder couldn't replace the '|cas_version|' variable with a value
    inside a code block.
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    cellarium_cas_version = get_tag()

    with open(f"{dir_path}/modules/installation.rst", "r") as f:
        content = f.read()

    # Replace the placeholder with the version
    content = content.replace("|cas_version|", cellarium_cas_version)

    # Write the processed content to the actual .rst file
    with open(f"{dir_path}/modules/installation_.rst", "w") as f:
        f.write(content)


substitute_installation_version()
