# Sphinx configuration for SPECWIZARD.
# Build locally with:   pip install -r docs/requirements.txt && make -C docs html
# Full option list: https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from importlib.metadata import PackageNotFoundError, version as _pkg_version

# Make the package importable without installing it (autodoc imports it).
sys.path.insert(0, os.path.abspath("../.."))

# -- Project information -----------------------------------------------------
project = "SPECWIZARD"
author = "Andrés Aramburo-García and Tom Theuns"
copyright = "2024–2026, " + author

try:  # single source of truth: the version in pyproject.toml
    release = _pkg_version("specwizard")
except PackageNotFoundError:
    import tomllib
    with open(os.path.join(os.path.dirname(__file__), "../../pyproject.toml"), "rb") as f:
        release = tomllib.load(f)["project"]["version"]
version = ".".join(release.split(".")[:2])

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",       # Google-style docstrings (what the code uses)
    "sphinx.ext.viewcode",       # "[source]" links in the API pages
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx_design",             # tabs, cards, grids (replaces sphinx_tabs)
    "sphinx_copybutton",
    "myst_nb",                   # notebooks + Markdown (replaces nbsphinx)
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "**.ipynb_checkpoints", "tutorials/configs"]
source_suffix = {".rst": "restructuredtext", ".md": "myst-nb", ".ipynb": "myst-nb"}

# -- autodoc / autosummary ---------------------------------------------------
autosummary_generate = True
autodoc_default_options = {"members": True, "show-inheritance": True}
autodoc_member_order = "bysource"
autoclass_content = "both"          # class docstring + __init__ docstring
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_use_ivar = True             # avoid duplicate attribute entries

# Heavy / simulation-specific dependencies are mocked so the API reference can
# be built on Read the Docs without installing them.
autodoc_mock_imports = [
    "swiftsimio", "pyread_eagle", "hydrangea", "mendeleev", "roman",
    "bs4", "requests", "html5lib", "astropy",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "unyt": ("https://unyt.readthedocs.io/en/stable/", None),
}

# -- Notebooks (myst-nb) -----------------------------------------------------
# The tutorials need simulation data that only lives on COSMA, so they cannot be
# executed on Read the Docs. Run them locally, commit them *with* outputs, and
# Sphinx renders the stored outputs.
nb_execution_mode = "off"
nb_merge_streams = True
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence", "deflist"]
myst_heading_anchors = 3

# -- HTML output -------------------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = f"SPECWIZARD {release}"
html_theme_options = {"navigation_depth": 3, "collapse_navigation": False}
