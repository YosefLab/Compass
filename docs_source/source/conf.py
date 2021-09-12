# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

import sys
from pathlib import Path

HERE = Path(__file__).parent
sys.path[:0] = [str(HERE.parent), str(HERE / "extensions")]


# -- Project information -----------------------------------------------------

project = 'Compass'
copyright = '2021, Yosef Lab'
author = 'Yosef Lab'

# The full version, including alpha/beta/rc tags
release = '2021'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'nbsphinx',
    'nbsphinx_link'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'pydata_sphinx_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


nbsphinx_prolog = r"""
.. raw:: html

{{% set docname = env.doc2path(env.docname, base=None).split("/")[-1].replace("nblink", "ipynb") %}}

.. raw:: html

    <style>
        p {{
            margin-bottom: 0.5rem;
        }}
    </style>

.. raw:: html

    <div class="admonition note">
        <p class="admonition-title">Note</p>
        <p>
        This page was generated from
        <a class="reference external" href="https://github.com/YosefLab/Compass/blob/docs/{docname}">{docname}</a>.
         Interactive online version:
        <span style="white-space: nowrap;"><a href="https://colab.research.google.com/github/YosefLab/Compass/blob/docs/{docname}"><img alt="Colab badge" src="https://colab.research.google.com/assets/colab-badge.svg" style="vertical-align:text-bottom"></a>.</span>
        </p>
    </div>
""".format(docname="{{ docname|e }}"
)

#Adds blurb to open in google collab automatically. Need to add way to add dataset to google collab noteboo though.
#Interactive online version:
#<span style="white-space: nowrap;"><a href="https://colab.research.google.com/github/YosefLab/Compass/blob/docs/{docname}"><img alt="Colab badge" src="https://colab.research.google.com/assets/colab-badge.svg" style="vertical-align:text-bottom"></a>.</span>
#