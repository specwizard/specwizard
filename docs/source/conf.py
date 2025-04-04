# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
sys.path.insert(0,os.path.abspath('../..'))

project = 'specwizard'
copyright = '2024, Andres Aramburo-Garcia and Tom Theuns'
author = 'Andres Aramburo and Tom Theuns'
release = '1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc','sphinx.ext.napoleon','sphinx_tabs.tabs','nbsphinx','matplotlib.sphinxext.mathmpl',
          'matplotlib.sphinxext.plot_directive',
          'sphinx.ext.doctest','sphinx_math_dollar','sphinx.ext.mathjax','myst_nb'
        ]
# Napoleon settings
napoleon_google_docstring = True

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
