# docs/conf.py — Sphinx configuration for LACAN
#
# Build HTML docs:
#   cd docs && make html
# Or from the repo root:
#   sphinx-build -b html docs docs/_build/html
#
# Requirements:
#   pip install sphinx sphinx-rtd-theme

import os
import sys

# Make the lacan package importable without installing it
sys.path.insert(0, os.path.abspath('.'))

# ---------------------------------------------------------------------------
# Project metadata
# ---------------------------------------------------------------------------

project = "LACAN"
copyright = "2026, Wim Dehaen"
author = "Wim Dehaen"
release = "1.0"

# ---------------------------------------------------------------------------
# Extensions
# ---------------------------------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",       # pull docstrings from source
    "sphinx.ext.napoleon",      # parse NumPy / Google style docstrings
    "sphinx.ext.viewcode",      # add [source] links to rendered pages
    "sphinx.ext.intersphinx",   # cross-link to Python / RDKit docs
    "sphinx.ext.autosummary",   # auto-generate summary tables
]

# autodoc: show both class docstring and __init__ docstring
autoclass_content = "both"

# autodoc: document members in source order, not alphabetical
autodoc_member_order = "bysource"

# autodoc: show type hints in the signature (not repeated in body)
autodoc_typehints = "signature"

# napoleon: accept both NumPy-style and Google-style docstrings
napoleon_numpy_docstring = True
napoleon_google_docstring = True
napoleon_use_param = True
napoleon_use_rtype = True

# autosummary: automatically generate stub .rst files
autosummary_generate = True

# intersphinx: cross-reference to external projects
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

# ---------------------------------------------------------------------------
# HTML output
# ---------------------------------------------------------------------------

html_theme = "sphinx_rtd_theme"

html_theme_options = {
    "navigation_depth": 4,
    "titles_only": False,
}

# Paths relative to docs/
html_static_path = ["_static"]
templates_path = ["_templates"]

# Files to ignore when looking for source files
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# ---------------------------------------------------------------------------
# Source files
# ---------------------------------------------------------------------------

# The master document (entry point)
master_doc = "index"
